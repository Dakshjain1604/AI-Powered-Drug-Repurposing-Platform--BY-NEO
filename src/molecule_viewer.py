"""
Molecular Visualization Module
Fetches molecular structures from PubChem and renders with stmol
"""

import requests
import logging
from typing import Optional, Dict
import streamlit as st

try:
    import stmol
    from stmol import showmol
    import py3Dmol
    STMOL_AVAILABLE = True
except ImportError:
    STMOL_AVAILABLE = False
    logging.warning("stmol not available, molecular visualization disabled")

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class MoleculeViewer:
    """
    Fetch and visualize molecular structures
    """
    
    def __init__(self, pubchem_base_url: str = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"):
        """
        Initialize molecule viewer
        
        Args:
            pubchem_base_url: Base URL for PubChem API
        """
        self.pubchem_base_url = pubchem_base_url
        self.timeout = 10
    
    def get_molecule_data(self, drug_name: str) -> Optional[Dict]:
        """
        Get molecular data for a drug from PubChem
        
        Args:
            drug_name: Name of the drug
            
        Returns:
            Dictionary with molecular data or None if not found
        """
        logger.info(f"Fetching molecular data for {drug_name}")
        
        # Get compound ID
        cid = self._get_compound_id(drug_name)
        if not cid:
            logger.warning(f"Could not find compound ID for {drug_name}")
            return None
        
        # Get molecular properties
        properties = self._get_compound_properties(cid)
        
        # Get SMILES
        smiles = self._get_smiles(cid)
        
        # Get 3D structure (SDF format)
        sdf = self._get_sdf(cid)
        
        return {
            "drug_name": drug_name,
            "cid": cid,
            "smiles": smiles,
            "sdf": sdf,
            "properties": properties,
            "pubchem_url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
        }
    
    def _get_compound_id(self, drug_name: str) -> Optional[str]:
        """
        Get PubChem compound ID (CID) from drug name
        """
        try:
            url = f"{self.pubchem_base_url}/compound/name/{drug_name}/cids/JSON"
            response = requests.get(url, timeout=self.timeout)
            
            if response.status_code == 200:
                data = response.json()
                cids = data.get("IdentifierList", {}).get("CID", [])
                if cids:
                    return str(cids[0])
            
            return None
            
        except Exception as e:
            logger.error(f"Error getting compound ID: {e}")
            return None
    
    def _get_compound_properties(self, cid: str) -> Dict:
        """
        Get compound properties from PubChem
        """
        try:
            url = f"{self.pubchem_base_url}/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,CanonicalSMILES,IUPACName/JSON"
            response = requests.get(url, timeout=self.timeout)
            
            if response.status_code == 200:
                data = response.json()
                properties = data.get("PropertyTable", {}).get("Properties", [{}])[0]
                return {
                    "molecular_formula": properties.get("MolecularFormula", "N/A"),
                    "molecular_weight": properties.get("MolecularWeight", "N/A"),
                    "iupac_name": properties.get("IUPACName", "N/A")
                }
            
            return {}
            
        except Exception as e:
            logger.error(f"Error getting compound properties: {e}")
            return {}
    
    def _get_smiles(self, cid: str) -> Optional[str]:
        """
        Get SMILES string for compound
        """
        try:
            url = f"{self.pubchem_base_url}/compound/cid/{cid}/property/CanonicalSMILES/JSON"
            response = requests.get(url, timeout=self.timeout)
            
            if response.status_code == 200:
                data = response.json()
                properties = data.get("PropertyTable", {}).get("Properties", [{}])[0]
                return properties.get("CanonicalSMILES")
            
            return None
            
        except Exception as e:
            logger.error(f"Error getting SMILES: {e}")
            return None
    
    def _get_sdf(self, cid: str) -> Optional[str]:
        """
        Get SDF (3D structure) for compound
        """
        try:
            url = f"{self.pubchem_base_url}/compound/cid/{cid}/record/SDF"
            response = requests.get(url, timeout=self.timeout)
            
            if response.status_code == 200:
                return response.text
            
            return None
            
        except Exception as e:
            logger.error(f"Error getting SDF: {e}")
            return None
    
    def render_molecule_streamlit(
        self, 
        molecule_data: Dict,
        width: int = 800,
        height: int = 600,
        style: str = "stick"
    ):
        """
        Render 3D molecule in Streamlit using stmol
        
        Args:
            molecule_data: Dictionary with molecular data (from get_molecule_data)
            width: Viewer width
            height: Viewer height
            style: Visualization style ('stick', 'sphere', 'cartoon')
        """
        if not STMOL_AVAILABLE:
            st.warning("Molecular visualization not available (stmol not installed)")
            return
        
        sdf = molecule_data.get("sdf")
        if not sdf:
            st.warning(f"No 3D structure available for {molecule_data.get('drug_name')}")
            return
        
        try:
            # Create 3D viewer
            viewer = py3Dmol.view(width=width, height=height)
            viewer.addModel(sdf, 'sdf')
            viewer.setStyle({style: {}})
            viewer.zoomTo()
            
            # Display in Streamlit
            showmol(viewer, height=height, width=width)
            
            logger.info(f"Rendered molecule: {molecule_data.get('drug_name')}")
            
        except Exception as e:
            logger.error(f"Error rendering molecule: {e}")
            st.error(f"Error rendering 3D structure: {e}")
    
    def display_molecule_info(self, molecule_data: Dict):
        """
        Display molecular information in Streamlit
        
        Args:
            molecule_data: Dictionary with molecular data
        """
        st.subheader(f"Molecular Information: {molecule_data.get('drug_name')}")
        
        properties = molecule_data.get("properties", {})
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.write("**Molecular Formula:**", properties.get("molecular_formula", "N/A"))
            st.write("**Molecular Weight:**", f"{properties.get('molecular_weight', 'N/A')} g/mol")
        
        with col2:
            st.write("**PubChem CID:**", molecule_data.get("cid", "N/A"))
            st.write("**SMILES:**", molecule_data.get("smiles", "N/A")[:50] + "..." if molecule_data.get("smiles") else "N/A")
        
        if molecule_data.get("pubchem_url"):
            st.write("**PubChem Link:**", molecule_data["pubchem_url"])


def get_and_display_molecule(drug_name: str, width: int = 800, height: int = 600):
    """
    Convenience function to fetch and display molecule in Streamlit
    
    Args:
        drug_name: Name of the drug
        width: Viewer width
        height: Viewer height
    """
    viewer = MoleculeViewer()
    
    with st.spinner(f"Fetching molecular data for {drug_name}..."):
        molecule_data = viewer.get_molecule_data(drug_name)
    
    if molecule_data:
        viewer.display_molecule_info(molecule_data)
        st.markdown("---")
        st.subheader("3D Structure")
        viewer.render_molecule_streamlit(molecule_data, width, height)
    else:
        st.warning(f"Could not retrieve molecular data for {drug_name}")
