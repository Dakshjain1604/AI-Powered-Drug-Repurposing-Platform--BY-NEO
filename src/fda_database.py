"""
FDA Drug Database Module
Validates drug approval status using local lookup and OpenFDA API
"""

import requests
from typing import Optional, Dict
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Simplified local database of known FDA-approved drugs (subset for common drugs)
# In production, this would be comprehensive or fully API-based
FDA_APPROVED_DRUGS_LOCAL = {
    "pirfenidone": {"approval_year": 2014, "indication": "Idiopathic pulmonary fibrosis"},
    "nintedanib": {"approval_year": 2014, "indication": "Idiopathic pulmonary fibrosis"},
    "sildenafil": {"approval_year": 1998, "indication": "Erectile dysfunction, Pulmonary arterial hypertension"},
    "losartan": {"approval_year": 1995, "indication": "Hypertension, diabetic nephropathy"},
    "metformin": {"approval_year": 1994, "indication": "Type 2 diabetes"},
    "simvastatin": {"approval_year": 1991, "indication": "Hypercholesterolemia"},
    "atorvastatin": {"approval_year": 1996, "indication": "Hypercholesterolemia"},
    "aspirin": {"approval_year": 1950, "indication": "Pain, inflammation, cardiovascular disease prevention"},
    "ibuprofen": {"approval_year": 1974, "indication": "Pain, inflammation, fever"},
    "prednisone": {"approval_year": 1955, "indication": "Inflammatory and autoimmune conditions"},
    "azathioprine": {"approval_year": 1968, "indication": "Immunosuppression, autoimmune diseases"},
    "mycophenolate": {"approval_year": 1995, "indication": "Immunosuppression, transplant rejection prevention"},
    "colchicine": {"approval_year": 1961, "indication": "Gout, familial Mediterranean fever"},
    "pentoxifylline": {"approval_year": 1984, "indication": "Peripheral vascular disease"},
    "n-acetylcysteine": {"approval_year": 1963, "indication": "Acetaminophen overdose, mucolytic"},
    "vitamin e": {"approval_year": 1989, "indication": "Vitamin E deficiency"},
    "ursodiol": {"approval_year": 1997, "indication": "Primary biliary cholangitis, gallstone dissolution"},
    "obeticholic acid": {"approval_year": 2016, "indication": "Primary biliary cholangitis"},
    "spironolactone": {"approval_year": 1960, "indication": "Heart failure, hypertension, edema"},
    "furosemide": {"approval_year": 1966, "indication": "Edema, hypertension"},
    "propranolol": {"approval_year": 1967, "indication": "Hypertension, angina, arrhythmia"},
    "methotrexate": {"approval_year": 1953, "indication": "Cancer, rheumatoid arthritis, psoriasis"},
    "adalimumab": {"approval_year": 2002, "indication": "Rheumatoid arthritis, Crohn's disease"},
    "infliximab": {"approval_year": 1998, "indication": "Crohn's disease, rheumatoid arthritis"},
    "tocilizumab": {"approval_year": 2010, "indication": "Rheumatoid arthritis"},
    "rituximab": {"approval_year": 1997, "indication": "Non-Hodgkin lymphoma, rheumatoid arthritis"},
    "cyclosporine": {"approval_year": 1983, "indication": "Immunosuppression, transplant rejection prevention"},
    "tacrolimus": {"approval_year": 1994, "indication": "Immunosuppression, transplant rejection prevention"},
    "sirolimus": {"approval_year": 1999, "indication": "Immunosuppression, transplant rejection prevention"},
    "everolimus": {"approval_year": 2009, "indication": "Cancer, transplant rejection prevention"},
    "sorafenib": {"approval_year": 2005, "indication": "Hepatocellular carcinoma, renal cell carcinoma"},
    "regorafenib": {"approval_year": 2012, "indication": "Colorectal cancer, hepatocellular carcinoma"},
    "vitamin d": {"approval_year": 1990, "indication": "Vitamin D deficiency, osteoporosis prevention"},
}

class FDADatabase:
    """
    FDA drug approval validation using local lookup and OpenFDA API
    """
    
    def __init__(self, use_api: bool = True, api_timeout: int = 10):
        """
        Initialize FDA database
        
        Args:
            use_api: Whether to use OpenFDA API for validation
            api_timeout: Timeout for API requests in seconds
        """
        self.use_api = use_api
        self.api_timeout = api_timeout
        self.openfda_base_url = "https://api.fda.gov/drug/drugsfda.json"
    
    def is_fda_approved(self, drug_name: str) -> bool:
        """
        Check if a drug is FDA-approved
        
        Args:
            drug_name: Name of the drug
            
        Returns:
            True if FDA-approved, False otherwise
        """
        drug_name_lower = drug_name.lower().strip()
        
        # Check local database first
        if drug_name_lower in FDA_APPROVED_DRUGS_LOCAL:
            return True
        
        # Try API if enabled
        if self.use_api:
            try:
                return self._check_openfda_api(drug_name)
            except Exception as e:
                logger.warning(f"OpenFDA API check failed for {drug_name}: {e}")
        
        return False
    
    def get_drug_info(self, drug_name: str) -> Optional[Dict]:
        """
        Get detailed information about an FDA-approved drug
        
        Args:
            drug_name: Name of the drug
            
        Returns:
            Dictionary with drug information or None if not found
        """
        drug_name_lower = drug_name.lower().strip()
        
        # Check local database
        if drug_name_lower in FDA_APPROVED_DRUGS_LOCAL:
            info = FDA_APPROVED_DRUGS_LOCAL[drug_name_lower].copy()
            info["drug_name"] = drug_name
            info["source"] = "local_database"
            return info
        
        # Try API if enabled
        if self.use_api:
            try:
                return self._get_openfda_info(drug_name)
            except Exception as e:
                logger.warning(f"OpenFDA API info retrieval failed for {drug_name}: {e}")
        
        return None
    
    def _check_openfda_api(self, drug_name: str) -> bool:
        """
        Check OpenFDA API for drug approval status
        
        Args:
            drug_name: Name of the drug
            
        Returns:
            True if found in OpenFDA, False otherwise
        """
        try:
            params = {
                "search": f'openfda.brand_name:"{drug_name}" OR openfda.generic_name:"{drug_name}"',
                "limit": 1
            }
            response = requests.get(self.openfda_base_url, params=params, timeout=self.api_timeout)
            
            if response.status_code == 200:
                data = response.json()
                return len(data.get("results", [])) > 0
            elif response.status_code == 404:
                return False
            else:
                logger.warning(f"OpenFDA API returned status {response.status_code}")
                return False
                
        except requests.exceptions.RequestException as e:
            logger.error(f"OpenFDA API request failed: {e}")
            return False
    
    def _get_openfda_info(self, drug_name: str) -> Optional[Dict]:
        """
        Get drug information from OpenFDA API
        
        Args:
            drug_name: Name of the drug
            
        Returns:
            Dictionary with drug information or None
        """
        try:
            params = {
                "search": f'openfda.brand_name:"{drug_name}" OR openfda.generic_name:"{drug_name}"',
                "limit": 1
            }
            response = requests.get(self.openfda_base_url, params=params, timeout=self.api_timeout)
            
            if response.status_code == 200:
                data = response.json()
                results = data.get("results", [])
                
                if results:
                    result = results[0]
                    openfda = result.get("openfda", {})
                    
                    return {
                        "drug_name": drug_name,
                        "approval_year": None,  # Not always available in API
                        "indication": ", ".join(openfda.get("pharm_class_epc", ["Unknown"])),
                        "source": "openfda_api",
                        "brand_names": openfda.get("brand_name", []),
                        "generic_names": openfda.get("generic_name", [])
                    }
            
            return None
            
        except requests.exceptions.RequestException as e:
            logger.error(f"OpenFDA API info retrieval failed: {e}")
            return None
    
    def add_local_drug(self, drug_name: str, approval_year: int, indication: str):
        """
        Add a drug to the local database (for extending coverage)
        
        Args:
            drug_name: Name of the drug
            approval_year: Year of FDA approval
            indication: Primary indication
        """
        drug_name_lower = drug_name.lower().strip()
        FDA_APPROVED_DRUGS_LOCAL[drug_name_lower] = {
            "approval_year": approval_year,
            "indication": indication
        }
        logger.info(f"Added {drug_name} to local FDA database")


# Convenience function for simple checks
def is_fda_approved(drug_name: str) -> bool:
    """
    Quick check if a drug is FDA-approved
    
    Args:
        drug_name: Name of the drug
        
    Returns:
        True if FDA-approved, False otherwise
    """
    db = FDADatabase(use_api=False)  # Use local only for speed
    return db.is_fda_approved(drug_name)
