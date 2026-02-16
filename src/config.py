"""
Configuration file for BioScript application
API keys, constants, and pathway definitions
"""

import os

# Try to load environment variables from .env file if dotenv is available
try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass  # dotenv not installed, will use os.getenv defaults

# API Configuration
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY", "")
NCBI_EMAIL = os.getenv("NCBI_EMAIL", "user@example.com")
NCBI_API_KEY = os.getenv("NCBI_API_KEY", "")  # Optional but recommended for higher rate limits

# PubMed/Literature Retrieval Settings
PUBMED_MAX_RESULTS = 50
PUBMED_DEFAULT_TIME_RANGE = 5  # years
PUBMED_RETRIES = 3
PUBMED_RETRY_DELAY = 2  # seconds

# Scoring Configuration
SCORING_WEIGHTS = {
    "pathway_hit": 40,  # Direct pathway interaction mentioned in literature
    "mechanism_alignment": 30,  # Mechanism aligns with fibrosis biology
    "clinical_evidence": 20,  # Safety/efficacy in similar conditions
    "market_availability": 10  # Availability and cost considerations
}

# LLM Configuration
LLM_MODEL = "gpt-4o-mini"  # or "gpt-3.5-turbo" for faster/cheaper option
LLM_TEMPERATURE = 0.3  # Lower for more deterministic scientific outputs
LLM_MAX_TOKENS = 2000

# Experiment Design Templates
CELL_MODELS = {
    "LX-2": "Human hepatic stellate cell line - gold standard for liver fibrosis",
    "THP-1": "Human monocytic cell line - for inflammation studies",
    "Primary hepatocytes": "Primary human hepatocytes - most physiologically relevant"
}

ASSAY_TYPES = {
    "Western blot": ["α-SMA", "Collagen-I", "Collagen-III", "TGF-β1"],
    "ELISA": ["TGF-β1", "IL-6", "IL-1β", "TNF-α"],
    "qPCR": ["COL1A1", "ACTA2", "TGFB1", "TIMP1", "MMP2"],
    "Immunofluorescence": ["α-SMA", "Collagen-I", "Vimentin"]
}

SUCCESS_METRICS = {
    "fibrotic_markers": "≥30% reduction in α-SMA and Collagen-I",
    "cytokines": "≥25% reduction in pro-fibrotic cytokines (TGF-β, IL-6)",
    "statistical_significance": "p < 0.05"
}

# Molecular Visualization Settings
PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
MOLECULE_VIEWER_SIZE = (800, 600)

# Report Generation Settings
PDF_REPORT_SECTIONS = [
    "Executive Summary",
    "Disease Background",
    "Literature Search Methodology",
    "Candidate Identification",
    "Top 3 Candidates",
    "Structured Hypotheses",
    "Experimental Protocols",
    "Molecular Structures",
    "References"
]

# Application Settings
APP_TITLE = "BioScript - Drug Repurposing Platform"
APP_ICON = "🧬"
DEFAULT_DISEASE = "Liver Fibrosis"
AVAILABLE_DISEASES = [
    "Liver Fibrosis",
    "Pulmonary Fibrosis",
    "Kidney Fibrosis",
    "Cardiac Fibrosis",
    "Systemic Sclerosis"
]

# Cache Settings
ENABLE_CACHING = True
CACHE_TTL = 3600  # seconds (1 hour)

# Logging Configuration
LOG_LEVEL = "INFO"
LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
