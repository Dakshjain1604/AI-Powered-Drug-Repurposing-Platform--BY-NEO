"""
Fibrotic Pathway Database
Curated dictionary of fibrotic targets and their biological roles
"""

FIBROTIC_PATHWAYS = {
    "TGF-β signaling": {
        "description": "Transforming Growth Factor-beta signaling pathway - primary driver of fibrosis",
        "key_molecules": ["TGF-β1", "TGF-β2", "TGF-β3", "TGFBR1", "TGFBR2"],
        "biological_role": "Promotes myofibroblast activation, ECM production, and epithelial-mesenchymal transition",
        "therapeutic_relevance": "High - central pathway in liver fibrosis",
        "weight": 40
    },
    "myofibroblast activation": {
        "description": "Activation of hepatic stellate cells into collagen-producing myofibroblasts",
        "key_molecules": ["α-SMA", "PDGF", "CTGF", "Vimentin"],
        "biological_role": "Drives excessive ECM deposition and scar tissue formation",
        "therapeutic_relevance": "High - primary effector cells in fibrosis",
        "weight": 35
    },
    "SMAD signaling": {
        "description": "SMAD proteins mediate TGF-β intracellular signaling",
        "key_molecules": ["SMAD2", "SMAD3", "SMAD4", "SMAD7"],
        "biological_role": "Transcriptional regulation of pro-fibrotic genes",
        "therapeutic_relevance": "High - downstream of TGF-β",
        "weight": 35
    },
    "ECM remodeling": {
        "description": "Extracellular matrix production and degradation imbalance",
        "key_molecules": ["Collagen-I", "Collagen-III", "Fibronectin", "MMPs", "TIMPs"],
        "biological_role": "Excessive collagen deposition leads to organ dysfunction",
        "therapeutic_relevance": "High - hallmark of fibrosis",
        "weight": 30
    },
    "inflammation": {
        "description": "Chronic inflammatory response driving fibrotic progression",
        "key_molecules": ["IL-6", "IL-1β", "TNF-α", "CCL2", "macrophages"],
        "biological_role": "Recruitment and activation of inflammatory cells perpetuates injury",
        "therapeutic_relevance": "Moderate - contributes to disease progression",
        "weight": 25
    },
    "oxidative stress": {
        "description": "Reactive oxygen species-mediated cellular damage",
        "key_molecules": ["ROS", "NOX enzymes", "SOD", "Catalase"],
        "biological_role": "Triggers stellate cell activation and promotes inflammation",
        "therapeutic_relevance": "Moderate - indirect driver",
        "weight": 20
    },
    "Wnt signaling": {
        "description": "Wnt/β-catenin pathway activation in fibrosis",
        "key_molecules": ["Wnt3a", "β-catenin", "GSK-3β"],
        "biological_role": "Promotes myofibroblast proliferation and survival",
        "therapeutic_relevance": "Emerging target",
        "weight": 15
    },
    "hedgehog signaling": {
        "description": "Hedgehog pathway activation in hepatic stellate cells",
        "key_molecules": ["Shh", "Gli1", "Ptch1"],
        "biological_role": "Drives stellate cell activation and proliferation",
        "therapeutic_relevance": "Emerging target",
        "weight": 15
    }
}

FIBROTIC_MARKERS = {
    "α-SMA": {
        "full_name": "Alpha-smooth muscle actin",
        "marker_type": "myofibroblast activation",
        "detection_methods": ["Western blot", "Immunofluorescence", "qPCR"],
        "clinical_significance": "Gold standard marker for activated stellate cells"
    },
    "Collagen-I": {
        "full_name": "Type I Collagen",
        "marker_type": "ECM deposition",
        "detection_methods": ["Western blot", "Sirius Red staining", "qPCR", "ELISA"],
        "clinical_significance": "Primary structural component of fibrotic tissue"
    },
    "Collagen-III": {
        "full_name": "Type III Collagen",
        "marker_type": "ECM deposition",
        "detection_methods": ["Western blot", "qPCR"],
        "clinical_significance": "Early fibrosis marker"
    },
    "TGF-β1": {
        "full_name": "Transforming Growth Factor-beta 1",
        "marker_type": "pro-fibrotic cytokine",
        "detection_methods": ["ELISA", "Western blot", "qPCR"],
        "clinical_significance": "Master regulator of fibrosis"
    },
    "TIMP-1": {
        "full_name": "Tissue Inhibitor of Metalloproteinases 1",
        "marker_type": "ECM remodeling inhibitor",
        "detection_methods": ["ELISA", "Western blot"],
        "clinical_significance": "Prevents ECM degradation, promotes fibrosis"
    }
}

def get_pathway_weight(pathway_name: str) -> int:
    """
    Get the scoring weight for a given pathway
    
    Args:
        pathway_name: Name of the fibrotic pathway
        
    Returns:
        Weight value (0-40) for scoring
    """
    pathway = FIBROTIC_PATHWAYS.get(pathway_name)
    if pathway:
        return pathway.get("weight", 0)
    return 0

def get_all_pathway_names() -> list:
    """
    Get list of all pathway names
    
    Returns:
        List of pathway names
    """
    return list(FIBROTIC_PATHWAYS.keys())

def get_pathway_info(pathway_name: str) -> dict:
    """
    Get detailed information about a pathway
    
    Args:
        pathway_name: Name of the pathway
        
    Returns:
        Dictionary with pathway details
    """
    return FIBROTIC_PATHWAYS.get(pathway_name, {})

def search_pathways_by_molecule(molecule: str) -> list:
    """
    Find pathways associated with a specific molecule
    
    Args:
        molecule: Molecule name (e.g., 'TGF-β1', 'α-SMA')
        
    Returns:
        List of pathway names where molecule is involved
    """
    matching_pathways = []
    molecule_lower = molecule.lower()
    
    for pathway_name, pathway_data in FIBROTIC_PATHWAYS.items():
        key_molecules = [m.lower() for m in pathway_data.get("key_molecules", [])]
        if any(molecule_lower in km or km in molecule_lower for km in key_molecules):
            matching_pathways.append(pathway_name)
    
    return matching_pathways
