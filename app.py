"""
BioResearch Agent Workbench - Full-Stack Research Dashboard
Data-dense interface for transparent biomedical research workflows
"""

import streamlit as st
import pandas as pd
import json
import logging
from datetime import datetime
import sys
import os
import plotly.graph_objects as go
import plotly.express as px
from io import BytesIO

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from src.config import (
    APP_TITLE, APP_ICON, DEFAULT_DISEASE, AVAILABLE_DISEASES,
    PUBMED_MAX_RESULTS, PUBMED_DEFAULT_TIME_RANGE, OPENAI_API_KEY, NCBI_EMAIL
)
from src.pubmed_fetcher import PubMedFetcher
from src.scientist_agent import ScientistAgent
from src.molecule_viewer import MoleculeViewer, get_and_display_molecule
from src.report_generator import ReportGenerator

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

# Page configuration
st.set_page_config(
    page_title="BioResearch Agent Workbench",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Dashboard-style CSS with muted professional theme
st.markdown("""
<style>
    /* High-density dashboard layout */
    .block-container {
        max-width: 100%;
        padding-left: 1.5rem;
        padding-right: 1.5rem;
        padding-top: 1rem;
    }
    
    /* Enhanced metric cards for dashboard */
    [data-testid="stMetricValue"] {
        font-size: 2rem;
        font-weight: 700;
    }
    
    [data-testid="stMetricLabel"] {
        font-size: 0.85rem;
        font-weight: 600;
        text-transform: uppercase;
        letter-spacing: 0.5px;
        opacity: 0.75;
    }
    
    /* Professional button styling */
    .stButton > button {
        border-radius: 4px;
        font-weight: 600;
        transition: all 0.2s ease;
        padding: 0.6rem 1.4rem;
    }
    
    /* Tab styling for dense dashboard */
    .stTabs [data-baseweb="tab-list"] {
        gap: 4px;
    }
    
    .stTabs [data-baseweb="tab"] {
        height: 52px;
        border-radius: 4px 4px 0 0;
        padding: 0 24px;
        font-weight: 600;
        font-size: 0.95rem;
    }
    
    /* DataFrames optimized for readability */
    .dataframe {
        font-size: 0.85rem;
    }
    
    /* Compact expanders */
    .streamlit-expanderHeader {
        border-radius: 4px;
        font-weight: 600;
        padding: 0.7rem 1rem;
        font-size: 0.9rem;
    }
    
    /* Clean layout */
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    
    /* Dashboard header styling */
    .dashboard-header {
        background: linear-gradient(135deg, rgba(90, 108, 125, 0.1), rgba(158, 174, 189, 0.1));
        padding: 1.5rem;
        border-radius: 8px;
        margin-bottom: 1.5rem;
    }
    
    /* Status badges */
    .status-badge {
        display: inline-block;
        padding: 0.3rem 0.8rem;
        border-radius: 4px;
        font-size: 0.85rem;
        font-weight: 600;
        margin: 0.2rem;
    }
    
    /* Reasoning trace terminal styling */
    .reasoning-trace {
        background-color: #1e1e1e;
        color: #d4d4d4;
        font-family: 'Monaco', 'Menlo', 'Ubuntu Mono', monospace;
        padding: 1rem;
        border-radius: 4px;
        font-size: 0.85rem;
        max-height: 500px;
        overflow-y: auto;
        line-height: 1.5;
    }
    
    .reasoning-step {
        margin-bottom: 0.5rem;
        padding: 0.3rem 0;
        border-left: 3px solid transparent;
        padding-left: 0.5rem;
    }
    
    .reasoning-step.initialization { border-left-color: #4CAF50; }
    .reasoning-step.extraction_start { border-left-color: #2196F3; }
    .reasoning-step.llm_query { border-left-color: #9C27B0; }
    .reasoning-step.fda_check { border-left-color: #FF9800; }
    .reasoning-step.scoring_candidate { border-left-color: #00BCD4; }
    .reasoning-step.error { border-left-color: #F44336; }
    
    .reasoning-timestamp {
        color: #858585;
        font-size: 0.75rem;
        margin-right: 0.5rem;
    }
    
    .reasoning-message {
        color: #d4d4d4;
    }
</style>

""", unsafe_allow_html=True)

# Initialize session state
def init_session_state():
    """Initialize all session state variables for research workbench"""
    defaults = {
        'abstracts': [],
        'all_extracted_drugs': [],
        'filtered_out_drugs': [],
        'extraction_logs': [],
        'fda_approved_candidates': [],
        'scored_candidates': [],
        'scoring_logs': [],
        'hypotheses': [],
        'protocols': [],
        'analysis_complete': False,
        'disease_focus': DEFAULT_DISEASE,
        'time_range': PUBMED_DEFAULT_TIME_RANGE,
        'repurposing_goal': "",
        'reasoning_trace': [],  # NEW: Store reasoning steps
        'reasoning_container': None,  # NEW: Container for live updates
    }
    
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value

init_session_state()


def main():
    """Main Research Workbench application"""
    
    # Dashboard header
    st.markdown('<div class="dashboard-header">', unsafe_allow_html=True)
    col1, col2 = st.columns([3, 1])
    with col1:
        st.title("🔬 BioResearch Agent Workbench")
        st.caption("Data-Driven Biomedical Research Platform • Transparent Pipeline • Full Data Access")
    with col2:
        if st.session_state.analysis_complete:
            st.markdown('<span class="status-badge" style="background-color: rgba(76, 175, 80, 0.2); color: #2e7d32;">✓ Analysis Complete</span>', unsafe_allow_html=True)
        else:
            st.markdown('<span class="status-badge" style="background-color: rgba(33, 150, 243, 0.2); color: #1565c0;">● Ready</span>', unsafe_allow_html=True)
    st.markdown('</div>', unsafe_allow_html=True)
    
    # Collapsible instructions
    with st.expander("📖 Research Workbench Guide", expanded=False):
        st.markdown("""
        **Research Workflow:**
        1. **Configure**: Set disease focus, time range, and research objective in sidebar
        2. **Execute**: Provide API credentials and run full analysis pipeline
        3. **Inspect**: View raw data in Data Inspector tab (all abstracts, extraction logs)
        4. **Analyze**: Examine filtering/scoring logic in Analysis Engine tab
        5. **Synthesize**: Review top candidates, hypotheses, and protocols in Synthesis tab
        6. **Export**: Generate comprehensive PDF report
        
        **Key Features:**
        - 🔍 **Full Data Transparency**: See ALL intermediate steps (raw abstracts, filtered candidates, scoring breakdown)
        - 📊 **Interactive Visualizations**: Explore data with charts and tables
        - 🧬 **3D Molecular Structures**: View candidate drug structures
        - 📄 **Professional Reports**: Export findings to PDF
        """)
    
    st.divider()
    
    # Sidebar Configuration
    render_sidebar()
    
    # Main Content Area
    if st.session_state.analysis_complete:
        render_research_workbench()
    else:
        render_ready_state()


def render_sidebar():
    """Render research configuration sidebar"""
    with st.sidebar:
        st.header("⚙️ Research Configuration")
        
        # Primary inputs
        disease_focus = st.selectbox(
            "Disease Focus",
            options=AVAILABLE_DISEASES,
            index=AVAILABLE_DISEASES.index(st.session_state.disease_focus),
            help="Target disease for research analysis"
        )
        st.session_state.disease_focus = disease_focus
        
        time_range = st.slider(
            "Publication Window (years)",
            min_value=1,
            max_value=10,
            value=st.session_state.time_range,
            help="Analyze publications from the last N years"
        )
        st.session_state.time_range = time_range
        
        repurposing_goal = st.text_area(
            "Research Objective",
            value=st.session_state.repurposing_goal or f"Identify FDA-approved drugs that modulate fibrotic pathways and reduce {disease_focus.lower()} progression",
            height=90,
            help="Define specific research objective"
        )
        st.session_state.repurposing_goal = repurposing_goal
        
        st.divider()
        
        # Advanced settings
        with st.expander("🔧 Advanced Settings", expanded=False):
            max_results = st.number_input(
                "Max Abstracts",
                min_value=10,
                max_value=200,
                value=PUBMED_MAX_RESULTS,
                step=10
            )
            
            min_score = st.slider(
                "Min Candidate Score",
                min_value=0,
                max_value=100,
                value=30
            )
        
        st.divider()
        
        # API Configuration
        st.header("🔑 API Configuration")
        
        openai_key_env = os.getenv("OPENAI_API_KEY", "")
        ncbi_email_env = os.getenv("NCBI_EMAIL", "")
        
        if openai_key_env:
            st.success("✓ OpenAI API Key (env)")
            openai_key = openai_key_env
            st.text_input(
                "OpenAI API Key",
                value="••••" + openai_key_env[-4:] if len(openai_key_env) > 4 else "••••",
                disabled=True
            )
        else:
            st.info("OpenAI API Key required")
            openai_key = st.text_input(
                "OpenAI API Key",
                type="password",
                placeholder="sk-..."
            )
        
        if ncbi_email_env and ncbi_email_env != "user@example.com":
            st.success("✓ NCBI Email (env)")
            ncbi_email = ncbi_email_env
            st.text_input(
                "NCBI Email",
                value=ncbi_email_env,
                disabled=True
            )
        else:
            st.info("NCBI Email required")
            ncbi_email = st.text_input(
                "NCBI Email",
                placeholder="researcher@institution.edu"
            )
        
        st.divider()
        
        # Action buttons
        col1, col2 = st.columns([3, 1])
        with col1:
            run_analysis = st.button(
                "▶ Start Analysis",
                type="primary",
                use_container_width=True,
                disabled=st.session_state.analysis_complete
            )
        with col2:
            if st.session_state.analysis_complete:
                if st.button("↻", help="Reset", use_container_width=True):
                    for key in list(st.session_state.keys()):
                        if key not in ['disease_focus', 'time_range', 'repurposing_goal']:
                            del st.session_state[key]
                    init_session_state()
                    st.rerun()
        
        if run_analysis:
            if not openai_key or (not openai_key_env and len(openai_key) < 20):
                st.error("❌ Valid OpenAI API Key required")
            elif not ncbi_email or '@' not in ncbi_email:
                st.error("❌ Valid NCBI Email required")
            else:
                run_full_analysis(
                    disease_focus, time_range, repurposing_goal,
                    max_results, openai_key, ncbi_email, min_score
                )


def render_ready_state():
    """Render ready state with preview metrics"""
    st.info("🔬 **Research Workbench Ready** • Configure parameters in sidebar and start analysis to populate the dashboard")
    
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Status", "Ready", help="System ready for analysis")
    
    with col2:
        st.metric("Disease Focus", st.session_state.disease_focus)
    
    with col3:
        st.metric("Time Window", f"{st.session_state.time_range} yrs")
    
    with col4:
        st.metric("Pipeline Stages", "5", help="Literature → Extract → Score → Hypothesize → Protocol")




def reasoning_callback(step_type: str, message: str, data: dict):
    """Callback function to capture reasoning steps in real-time"""
    if 'reasoning_trace' not in st.session_state:
        st.session_state.reasoning_trace = []
    
    trace_entry = {
        "timestamp": datetime.now().strftime("%H:%M:%S"),
        "step_type": step_type,
        "message": message,
        "data": data
    }
    st.session_state.reasoning_trace.append(trace_entry)
    
    # Update live container if it exists
    if st.session_state.get('reasoning_container'):
        try:
            with st.session_state.reasoning_container:
                render_reasoning_step(trace_entry)
        except:
            pass  # Container may not be ready


def render_reasoning_step(step: dict):
    """Render a single reasoning step"""
    timestamp = step.get("timestamp", "")
    message = step.get("message", "")
    step_type = step.get("step_type", "")
    
    st.markdown(
        f'<div class="reasoning-step {step_type}">'
        f'<span class="reasoning-timestamp">[{timestamp}]</span>'
        f'<span class="reasoning-message">{message}</span>'
        f'</div>',
        unsafe_allow_html=True
    )


def render_reasoning_trace_panel():
    """Render full reasoning trace panel"""
    st.subheader("🧠 Agent Reasoning Trace")
    st.caption("Real-time view of agent's internal decision-making process")
    
    if not st.session_state.reasoning_trace:
        st.info("No reasoning trace available. Run analysis to see agent thinking.")
        return
    
    # Container for reasoning trace
    with st.container():
        st.markdown('<div class="reasoning-trace">', unsafe_allow_html=True)
        
        for step in st.session_state.reasoning_trace:
            render_reasoning_step(step)
        
        st.markdown('</div>', unsafe_allow_html=True)
    
    # Expandable detailed view
    with st.expander("📊 View Detailed Trace Data (JSON)", expanded=False):
        st.json(st.session_state.reasoning_trace)



def run_full_analysis(
    disease_focus: str,
    time_range: int,
    repurposing_goal: str,
    max_results: int,
    openai_key: str,
    ncbi_email: str,
    min_score: int
):
    """Execute complete research pipeline with full data tracking"""
    
    # Reset reasoning trace
    st.session_state.reasoning_trace = []
    
    # Create live reasoning feed container
    st.subheader("🧠 Live Agent Reasoning Feed")
    st.caption("Watch the agent's thinking process in real-time")
    
    reasoning_container = st.container()
    st.session_state.reasoning_container = reasoning_container
    
    with reasoning_container:
        st.markdown('<div class="reasoning-trace">', unsafe_allow_html=True)
    
    with st.status("🔄 Running Research Pipeline...", expanded=True) as status:
        try:
            # Step 1: Literature Retrieval
            st.write("**Stage 1/5:** Fetching literature from PubMed...")
            reasoning_callback("literature_start", "🔍 Initializing PubMed fetcher", {})
            
            fetcher = PubMedFetcher(email=ncbi_email)
            query = f"{disease_focus} fibrosis drug therapy treatment"
            
            reasoning_callback("literature_query", f"📚 Searching PubMed: {query}", {"query": query, "max_results": max_results})
            abstracts = fetcher.fetch_literature(query, max_results, time_range)
            
            if not abstracts:
                reasoning_callback("error", "❌ No abstracts found", {})
                st.error("❌ No abstracts found. Adjust search parameters.")
                status.update(label="Analysis Failed", state="error")
                return
            
            st.session_state.abstracts = abstracts
            reasoning_callback("literature_complete", f"✅ Retrieved {len(abstracts)} abstracts from PubMed", {"count": len(abstracts)})
            st.write(f"✅ Retrieved **{len(abstracts)}** abstracts from PubMed")
            
            # Step 2: Candidate Extraction (with reasoning trace)
            st.write("**Stage 2/5:** Extracting drug candidates with full transparency...")
            reasoning_callback("extraction_pipeline_start", "🔬 Initializing candidate extraction pipeline", {})
            
            agent = ScientistAgent(
                openai_api_key=openai_key,
                reasoning_callback=reasoning_callback  # Pass callback to agent
            )
            
            extraction_result = agent.extract_candidates_from_text(abstracts, disease_focus)
            
            # Store ALL intermediate data
            st.session_state.fda_approved_candidates = extraction_result.get('fda_approved', [])
            st.session_state.all_extracted_drugs = extraction_result.get('all_extracted', [])
            st.session_state.filtered_out_drugs = extraction_result.get('filtered_out', [])
            st.session_state.extraction_logs = extraction_result.get('extraction_logs', [])
            
            if not st.session_state.fda_approved_candidates:
                reasoning_callback("warning", "⚠️ No FDA-approved candidates found", {})
                st.warning("⚠️ No FDA-approved candidates found in abstracts.")
                status.update(label="Analysis Incomplete", state="error")
                return
            
            reasoning_callback("extraction_summary", 
                f"✅ Extraction complete: {len(st.session_state.all_extracted_drugs)} total, {len(st.session_state.fda_approved_candidates)} FDA-approved",
                {
                    "total_extracted": len(st.session_state.all_extracted_drugs),
                    "fda_approved": len(st.session_state.fda_approved_candidates),
                    "filtered_out": len(st.session_state.filtered_out_drugs)
                }
            )
            st.write(f"✅ Extracted **{len(st.session_state.all_extracted_drugs)}** total drugs "
                    f"(FDA-approved: **{len(st.session_state.fda_approved_candidates)}**, "
                    f"Filtered: **{len(st.session_state.filtered_out_drugs)}**)")
            
            # Step 3: Scoring (with reasoning trace)
            st.write("**Stage 3/5:** Scoring candidates with pathway analysis...")
            reasoning_callback("scoring_pipeline_start", "📊 Starting candidate scoring pipeline", {})
            
            scoring_result = agent.score_candidates(
                st.session_state.fda_approved_candidates, 
                abstracts
            )
            
            st.session_state.scored_candidates = scoring_result.get('scored_candidates', [])
            st.session_state.scoring_logs = scoring_result.get('scoring_logs', [])
            
            # Apply score threshold
            st.session_state.scored_candidates = [
                c for c in st.session_state.scored_candidates 
                if c.get('score', 0) >= min_score
            ]
            
            if not st.session_state.scored_candidates:
                reasoning_callback("warning", f"⚠️ No candidates above threshold {min_score}", {"threshold": min_score})
                st.warning(f"⚠️ No candidates scored above {min_score}/100. Lower threshold.")
                status.update(label="Analysis Incomplete", state="error")
                return
            
            reasoning_callback("scoring_complete", 
                f"✅ Scored {len(st.session_state.scored_candidates)} candidates (threshold: {min_score})",
                {"scored_count": len(st.session_state.scored_candidates), "threshold": min_score}
            )
            st.write(f"✅ Scored **{len(st.session_state.scored_candidates)}** candidates (threshold: {min_score}/100)")
            
            # Step 4: Hypothesis Generation
            st.write("**Stage 4/5:** Generating structured hypotheses...")
            reasoning_callback("hypothesis_pipeline_start", "💡 Starting hypothesis generation", {})
            
            top_3 = st.session_state.scored_candidates[:3]
            hypotheses = agent.generate_structured_hypothesis(top_3, disease_focus)
            st.session_state.hypotheses = hypotheses
            
            reasoning_callback("hypothesis_complete", f"✅ Generated {len(hypotheses)} hypotheses", {"count": len(hypotheses)})
            st.write(f"✅ Generated **{len(hypotheses)}** evidence-based hypotheses")
            
            # Step 5: Protocol Design
            st.write("**Stage 5/5:** Designing experimental protocols...")
            reasoning_callback("protocol_pipeline_start", "🧪 Starting experimental protocol design", {})
            
            protocols = []
            for candidate in top_3:
                protocol = agent.design_experiment_protocol(candidate)
                protocols.append(protocol)
            
            st.session_state.protocols = protocols
            reasoning_callback("protocol_complete", f"✅ Designed {len(protocols)} protocols", {"count": len(protocols)})
            st.write(f"✅ Designed **{len(protocols)}** detailed experimental protocols")
            
            st.session_state.analysis_complete = True
            reasoning_callback("pipeline_complete", "🎉 Research pipeline completed successfully", {
                "total_abstracts": len(abstracts),
                "total_candidates": len(st.session_state.scored_candidates),
                "top_3": [c.get("drug_name", "Unknown") for c in top_3]
            })
            
            status.update(label="✅ Research Pipeline Complete", state="complete")
            
            with reasoning_container:
                st.markdown('</div>', unsafe_allow_html=True)
            
            st.success(f"🎉 Successfully analyzed **{len(abstracts)}** abstracts and identified **{len(st.session_state.scored_candidates)}** candidates")
            st.rerun()
            
        except Exception as e:
            reasoning_callback("error", f"❌ Pipeline error: {str(e)}", {"error": str(e)})
            st.error(f"❌ Pipeline Error: {str(e)}")
            logger.error(f"Analysis error: {e}", exc_info=True)
            status.update(label="Analysis Failed", state="error")
            
            with st.expander("🔍 Error Details"):
                st.code(str(e))

def render_research_workbench():
    """Render main research workbench with data-dense tabs"""
    
    # Dashboard metrics
    render_dashboard_metrics()
    
    st.divider()
    
    # Tab-based research workbench
    tabs = st.tabs([
        "📊 Data Inspector",
        "🔬 Analysis Engine",
        "🎯 Synthesis",
        "🧠 Reasoning Trace",
        "📄 Export & Visualization"
    ])
    
    with tabs[0]:
        render_data_inspector_tab()
    
    with tabs[1]:
        render_analysis_engine_tab()
    
    with tabs[2]:
        render_synthesis_tab()
    
    with tabs[3]:
        render_reasoning_trace_panel()
    
    with tabs[4]:
        render_export_tab()


def render_dashboard_metrics():
    """Render comprehensive dashboard metrics"""
    st.subheader("📈 Research Dashboard")
    
    col1, col2, col3, col4, col5, col6 = st.columns(6)
    
    abstracts_count = len(st.session_state.abstracts)
    total_extracted = len(st.session_state.all_extracted_drugs)
    fda_approved = len(st.session_state.fda_approved_candidates)
    filtered_out = len(st.session_state.filtered_out_drugs)
    scored_count = len(st.session_state.scored_candidates)
    top_score = st.session_state.scored_candidates[0].get('score', 0) if scored_count > 0 else 0
    
    with col1:
        st.metric("Abstracts", abstracts_count, help="Total PubMed abstracts retrieved")
    
    with col2:
        st.metric("Extracted", total_extracted, help="Total drugs extracted from abstracts")
    
    with col3:
        st.metric("FDA Approved", fda_approved, delta=f"-{filtered_out} filtered", help="FDA-approved candidates")
    
    with col4:
        st.metric("Scored", scored_count, help="Candidates meeting score threshold")
    
    with col5:
        st.metric("Top Score", f"{top_score:.0f}/100", help="Highest candidate score")
    
    with col6:
        top_drug = st.session_state.scored_candidates[0].get('drug_name', 'N/A') if scored_count > 0 else 'N/A'
        st.metric("Leader", top_drug[:15], help="Top-ranked candidate")


def render_data_inspector_tab():
    """Data Inspector: View ALL raw research data"""
    
    st.subheader("🔍 Data Inspector - Raw Research Data")
    st.caption("View all intermediate data for complete research transparency")
    
    # Sub-tabs for different data views
    data_tabs = st.tabs([
        "Literature Corpus",
        "Extraction Results",
        "Filtering Logs"
    ])
    
    # Tab 1: Literature Corpus
    with data_tabs[0]:
        st.markdown("#### 📚 Complete Literature Corpus")
        st.caption(f"Showing all {len(st.session_state.abstracts)} retrieved abstracts from PubMed")
        
        # Build comprehensive dataframe
        abstracts_data = []
        for i, abstract in enumerate(st.session_state.abstracts, 1):
            abstracts_data.append({
                '#': i,
                'PMID': abstract.get('pmid', 'N/A'),
                'Title': abstract.get('title', 'No title'),
                'Journal': abstract.get('journal', 'N/A'),
                'Date': abstract.get('publication_date', 'N/A'),
                'Authors': ', '.join(abstract.get('authors', [])[:3]) + ("..." if len(abstract.get('authors', [])) > 3 else ""),
                'Abstract': abstract.get('abstract_text', '')[:200] + "..."
            })
        
        df_abstracts = pd.DataFrame(abstracts_data)
        
        # Display with full width
        st.dataframe(
            df_abstracts,
            use_container_width=True,
            height=450,
            column_config={
                "#": st.column_config.NumberColumn(width="small"),
                "PMID": st.column_config.TextColumn(width="small"),
                "Title": st.column_config.TextColumn(width="large"),
                "Journal": st.column_config.TextColumn(width="medium"),
                "Date": st.column_config.TextColumn(width="small"),
                "Authors": st.column_config.TextColumn(width="medium"),
                "Abstract": st.column_config.TextColumn(width="large")
            }
        )
        
        # Download option
        csv_abstracts = pd.DataFrame([{
            'pmid': a.get('pmid'),
            'title': a.get('title'),
            'journal': a.get('journal'),
            'date': a.get('publication_date'),
            'doi': a.get('doi'),
            'authors': ', '.join(a.get('authors', [])),
            'abstract': a.get('abstract_text')
        } for a in st.session_state.abstracts]).to_csv(index=False)
        
        st.download_button(
            label="📥 Download Full Corpus (CSV)",
            data=csv_abstracts,
            file_name=f"literature_corpus_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
            mime="text/csv"
        )
    
    # Tab 2: Extraction Results
    with data_tabs[1]:
        st.markdown("#### 🧬 Drug Extraction Results")
        
        col1, col2 = st.columns([1, 1])
        
        with col1:
            st.metric("Total Extracted", len(st.session_state.all_extracted_drugs))
            st.caption("All drugs identified in abstracts (before FDA filtering)")
            
            if st.session_state.all_extracted_drugs:
                extraction_data = []
                for i, drug in enumerate(st.session_state.all_extracted_drugs, 1):
                    extraction_data.append({
                        '#': i,
                        'Drug Name': drug.get('drug_name', 'Unknown'),
                        'Mechanism': drug.get('mechanism', 'N/A')[:80] + "...",
                        'Pathways': ', '.join(drug.get('pathway_targets', [])[:2]),
                        'Evidence': drug.get('evidence_type', 'N/A'),
                        'Confidence': f"{drug.get('extraction_confidence', 0):.2f}"
                    })
                
                df_extracted = pd.DataFrame(extraction_data)
                st.dataframe(df_extracted, use_container_width=True, height=350)
        
        with col2:
            st.metric("FDA Approved", len(st.session_state.fda_approved_candidates))
            st.caption("Candidates that passed FDA approval verification")
            
            if st.session_state.fda_approved_candidates:
                fda_data = []
                for i, drug in enumerate(st.session_state.fda_approved_candidates, 1):
                    fda_data.append({
                        '#': i,
                        'Drug Name': drug.get('drug_name', 'Unknown'),
                        'Approval Year': drug.get('approval_year', 'N/A'),
                        'Original Indication': drug.get('indication', 'Unknown')[:50] + "...",
                        'Pathways': ', '.join(drug.get('pathway_targets', [])[:2])
                    })
                
                df_fda = pd.DataFrame(fda_data)
                st.dataframe(df_fda, use_container_width=True, height=350)
    
    # Tab 3: Filtering Logs
    with data_tabs[2]:
        st.markdown("#### 🔍 Filtering & Validation Logs")
        
        col1, col2 = st.columns([2, 1])
        
        with col1:
            st.metric("Filtered Out", len(st.session_state.filtered_out_drugs))
            st.caption("Drugs excluded due to FDA approval status")
            
            if st.session_state.filtered_out_drugs:
                filtered_data = []
                for i, drug in enumerate(st.session_state.filtered_out_drugs, 1):
                    filtered_data.append({
                        '#': i,
                        'Drug Name': drug.get('drug_name', 'Unknown'),
                        'Filter Reason': drug.get('filter_reason', 'Unknown'),
                        'Mechanism': drug.get('mechanism', 'N/A')[:60] + "...",
                        'Source PMID': drug.get('source_pmid', 'N/A')
                    })
                
                df_filtered = pd.DataFrame(filtered_data)
                st.dataframe(df_filtered, use_container_width=True, height=400)
        
        with col2:
            st.markdown("##### Extraction Timeline")
            if st.session_state.extraction_logs:
                for log in st.session_state.extraction_logs:
                    st.caption(f"**Timestamp:** {log.get('timestamp', 'N/A')[:19]}")
                    if 'total_extracted' in log:
                        st.caption(f"• Total Extracted: {log.get('total_extracted', 0)}")
                        st.caption(f"• Source Abstracts: {log.get('source_abstracts', 0)}")
                    if 'fda_approved_count' in log:
                        st.caption(f"• FDA Approved: {log.get('fda_approved_count', 0)}")
                        st.caption(f"• Filtered Out: {log.get('filtered_out_count', 0)}")
                    st.divider()


def render_analysis_engine_tab():
    """Analysis Engine: Visualize scoring and filtering logic"""
    
    st.subheader("🔬 Analysis Engine - Scoring & Ranking Logic")
    st.caption("Transparent view of candidate evaluation pipeline")
    
    # Funnel visualization
    st.markdown("#### 📉 Research Funnel")
    
    funnel_data = {
        'Stage': [
            'Abstracts Retrieved',
            'Drugs Extracted',
            'FDA Approved',
            'Passed Scoring',
            'Top Candidates'
        ],
        'Count': [
            len(st.session_state.abstracts),
            len(st.session_state.all_extracted_drugs),
            len(st.session_state.fda_approved_candidates),
            len(st.session_state.scored_candidates),
            min(3, len(st.session_state.scored_candidates))
        ]
    }
    
    fig_funnel = go.Figure(go.Funnel(
        y=funnel_data['Stage'],
        x=funnel_data['Count'],
        textinfo="value+percent initial",
        marker=dict(color=['#5a6c7d', '#6b7c8d', '#7c8d9d', '#8d9ead', '#9eaebd']),
        connector={"line": {"color": "#4a5a6a", "width": 2}}
    ))
    
    fig_funnel.update_layout(
        height=400,
        template="plotly_white",
        font=dict(family="sans-serif", size=11, color="#2c3e50"),
        paper_bgcolor='rgba(0,0,0,0)'
    )
    
    st.plotly_chart(fig_funnel, use_container_width=True)
    
    st.divider()
    
    # Scoring breakdown table
    st.markdown("#### 📊 Complete Scoring Breakdown")
    st.caption("Detailed scoring logic for all candidates")
    
    if st.session_state.scoring_logs:
        scoring_table_data = []
        for i, log in enumerate(st.session_state.scoring_logs, 1):
            breakdown = log.get('breakdown', {})
            scoring_table_data.append({
                'Rank': i,
                'Drug Name': log.get('drug_name', 'Unknown'),
                'Total Score': log.get('total_score', 0),
                'Pathway Hit': breakdown.get('pathway_hit', 0),
                'Mechanism': breakdown.get('mechanism_alignment', 0),
                'Clinical': breakdown.get('clinical_evidence', 0),
                'Availability': breakdown.get('market_availability', 0),
                'Target Pathways': ', '.join(log.get('pathway_targets', [])[:2])
            })
        
        df_scoring = pd.DataFrame(scoring_table_data)
        
        st.dataframe(
            df_scoring,
            use_container_width=True,
            height=450,
            column_config={
                "Rank": st.column_config.NumberColumn(width="small"),
                "Drug Name": st.column_config.TextColumn(width="medium"),
                "Total Score": st.column_config.ProgressColumn(
                    width="medium",
                    format="%d/100",
                    min_value=0,
                    max_value=100
                ),
                "Pathway Hit": st.column_config.ProgressColumn(width="small", format="%d", min_value=0, max_value=40),
                "Mechanism": st.column_config.ProgressColumn(width="small", format="%d", min_value=0, max_value=30),
                "Clinical": st.column_config.ProgressColumn(width="small", format="%d", min_value=0, max_value=20),
                "Availability": st.column_config.ProgressColumn(width="small", format="%d", min_value=0, max_value=10),
                "Target Pathways": st.column_config.TextColumn(width="large")
            }
        )
        
        # Download scoring data
        csv_scoring = df_scoring.to_csv(index=False)
        st.download_button(
            label="📥 Download Scoring Data (CSV)",
            data=csv_scoring,
            file_name=f"scoring_breakdown_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
            mime="text/csv"
        )
    
    st.divider()
    
    # Top candidates comparison
    st.markdown("#### 🏆 Top Candidates - Score Comparison")
    
    top_10 = st.session_state.scored_candidates[:10]
    
    muted_colors = ['#5a6c7d', '#6b7c8d', '#7c8d9d', '#8d9ead', '#9eaebd', 
                    '#afbecd', '#c0cfdd', '#d1dfed', '#e2effd', '#f3f9ff']
    
    fig_comparison = go.Figure()
    
    fig_comparison.add_trace(go.Bar(
        x=[c.get('drug_name', 'Unknown') for c in top_10],
        y=[c.get('score', 0) for c in top_10],
        marker=dict(
            color=muted_colors[:len(top_10)],
            line=dict(color='#4a5a6a', width=1)
        ),
        text=[f"{c.get('score', 0)}/100" for c in top_10],
        textposition='outside',
        hovertemplate='<b>%{x}</b><br>Score: %{y}/100<extra></extra>'
    ))
    
    fig_comparison.update_layout(
        xaxis_title="Candidate Drug",
        yaxis_title="Total Score",
        yaxis=dict(range=[0, 105]),
        height=450,
        template="plotly_white",
        hovermode='x',
        font=dict(family="sans-serif", size=11, color="#2c3e50"),
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)'
    )
    
    st.plotly_chart(fig_comparison, use_container_width=True)


def render_synthesis_tab():
    """Synthesis: Display top candidates, hypotheses, and protocols"""
    
    st.subheader("🎯 Synthesis - Research Outcomes")
    st.caption("Top-ranked candidates with structured hypotheses and experimental protocols")
    
    candidates = st.session_state.scored_candidates[:3]
    hypotheses = st.session_state.hypotheses
    protocols = st.session_state.protocols
    
    for i, (candidate, hypothesis, protocol) in enumerate(zip(candidates, hypotheses, protocols), 1):
        medal = "🥇" if i == 1 else "🥈" if i == 2 else "🥉"
        
        st.markdown(f"### {medal} Rank #{i} - {candidate.get('drug_name', 'Unknown')}")
        
        # Metrics row
        metric_col1, metric_col2, metric_col3, metric_col4, metric_col5 = st.columns(5)
        
        breakdown = candidate.get('score_breakdown', {})
        
        with metric_col1:
            st.metric("Total Score", f"{candidate.get('score', 0)}/100")
        
        with metric_col2:
            st.metric("Pathway", f"{breakdown.get('pathway_hit', 0)}/40")
        
        with metric_col3:
            st.metric("Mechanism", f"{breakdown.get('mechanism_alignment', 0)}/30")
        
        with metric_col4:
            st.metric("Clinical", f"{breakdown.get('clinical_evidence', 0)}/20")
        
        with metric_col5:
            st.metric("Availability", f"{breakdown.get('market_availability', 0)}/10")
        
        # Detailed tabs for each candidate
        cand_tabs = st.tabs(["📋 Overview", "💡 Hypothesis", "🧪 Protocol", "📊 Details",
        "🧠 Reasoning Trace"])
        
        with cand_tabs[0]:
            col1, col2 = st.columns([1, 1])
            
            with col1:
                st.markdown("**Drug Information**")
                st.write(f"• **FDA Approval Year:** {candidate.get('approval_year', 'Unknown')}")
                st.write(f"• **Original Indication:** {candidate.get('indication', 'Unknown')}")
                st.write(f"• **Extraction Confidence:** {candidate.get('extraction_confidence', 0):.1%}")
                st.write(f"• **Source PMID:** {candidate.get('source_pmid', 'N/A')}")
            
            with col2:
                st.markdown("**Target Pathways**")
                for pathway in candidate.get('pathway_targets', []):
                    st.write(f"• {pathway}")
                
                st.markdown("**Mechanism of Action**")
                st.caption(candidate.get('mechanism', 'No mechanism described')[:200] + "...")
        
        with cand_tabs[1]:
            st.info(f"**Hypothesis:** {hypothesis.get('hypothesis', 'No hypothesis generated')}")
            
            st.markdown("**Expected Biological Outcomes:**")
            for outcome in hypothesis.get('expected_outcomes', []):
                st.write(f"• {outcome}")
            
            st.markdown("**Biological Reasoning:**")
            st.write(hypothesis.get('reasoning', 'No reasoning provided'))
            
            st.metric("Confidence Score", f"{hypothesis.get('confidence_score', 0):.1%}")
        
        with cand_tabs[2]:
            col1, col2 = st.columns([1, 1])
            
            with col1:
                st.markdown("**Experimental Design**")
                st.write(f"• **Cell Model:** {protocol.get('model_system', 'N/A')}")
                st.write(f"• **Assay Type:** {protocol.get('assay_type', 'N/A')}")
                st.write(f"• **Drug Concentration:** {protocol.get('drug_concentration_range', 'N/A')}")
                st.write(f"• **Dosing Schedule:** {protocol.get('dosing_schedule', 'N/A')}")
                st.write(f"• **Timeline:** {protocol.get('assay_timeline', 'N/A')}")
            
            with col2:
                st.markdown("**Success Metrics**")
                for metric in protocol.get('success_metrics', []):
                    st.write(f"✓ {metric}")
                
                st.markdown("**Controls**")
                controls = protocol.get('controls', {})
                st.write(f"• **Positive:** {controls.get('positive', 'N/A')}")
                st.write(f"• **Negative:** {controls.get('negative', 'N/A')}")
                
                st.markdown("**Statistical Analysis**")
                st.caption(protocol.get('statistical_analysis', 'N/A'))
        
        with cand_tabs[3]:
            st.markdown("**Score Breakdown Analysis**")
            
            breakdown_df = pd.DataFrame({
                'Component': ['Pathway Hit', 'Mechanism Alignment', 'Clinical Evidence', 'Market Availability'],
                'Score': [
                    breakdown.get('pathway_hit', 0),
                    breakdown.get('mechanism_alignment', 0),
                    breakdown.get('clinical_evidence', 0),
                    breakdown.get('market_availability', 0)
                ],
                'Max Possible': [40, 30, 20, 10],
                'Percentage': [
                    f"{(breakdown.get('pathway_hit', 0)/40)*100:.1f}%",
                    f"{(breakdown.get('mechanism_alignment', 0)/30)*100:.1f}%",
                    f"{(breakdown.get('clinical_evidence', 0)/20)*100:.1f}%",
                    f"{(breakdown.get('market_availability', 0)/10)*100:.1f}%"
                ]
            })
            
            st.dataframe(breakdown_df, use_container_width=True, hide_index=True)
        
        st.divider()


def render_export_tab():
    """Export: Generate reports and visualize molecules"""
    
    st.subheader("📄 Export & Visualization")
    st.caption("Generate comprehensive research reports and view molecular structures")
    
    # Report generation
    st.markdown("#### 📋 Generate Research Report")
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        report_title = st.text_input(
            "Report Title",
            value=f"BioResearch Agent Report - {st.session_state.disease_focus}"
        )
        
        include_all = st.checkbox(
            "Include all candidates (not just top 3)",
            value=False
        )
        
        include_abstracts = st.checkbox(
            "Include abstract summaries",
            value=True
        )
    
    with col2:
        st.caption("**Report Sections:**")
        st.caption("✓ Executive Summary")
        st.caption("✓ Research Methodology")
        st.caption("✓ Data Analysis")
        st.caption("✓ Candidate Ranking")
        st.caption("✓ Hypotheses & Protocols")
        st.caption("✓ References")
    
    if st.button("📥 Generate PDF Report", type="primary"):
        with st.spinner("Generating comprehensive research report..."):
            try:
                generator = ReportGenerator()
                
                candidates_to_include = st.session_state.scored_candidates if include_all else st.session_state.scored_candidates[:3]
                abstracts_to_include = st.session_state.abstracts if include_abstracts else []
                
                os.makedirs('/root/aiCoScientist/data', exist_ok=True)
                
                filepath = f"/root/aiCoScientist/data/report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdf"
                
                success = generator.generate_report(
                    filepath=filepath,
                    disease_focus=st.session_state.disease_focus,
                    time_range=f"last {st.session_state.time_range} years",
                    repurposing_goal=st.session_state.repurposing_goal,
                    abstracts=abstracts_to_include,
                    candidates=candidates_to_include,
                    hypotheses=st.session_state.hypotheses,
                    protocols=st.session_state.protocols
                )
                
                if success:
                    st.success("✅ Report generated successfully!")
                    
                    with open(filepath, 'rb') as f:
                        pdf_data = f.read()
                    
                    st.download_button(
                        label="📥 Download PDF Report",
                        data=pdf_data,
                        file_name=f"bioresearch_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdf",
                        mime="application/pdf"
                    )
                    
                    st.info(f"Report includes: {len(candidates_to_include)} candidates • {len(st.session_state.hypotheses)} hypotheses • {len(st.session_state.protocols)} protocols")
                else:
                    st.error("❌ Error generating report")
                    
            except Exception as e:
                st.error(f"❌ Error: {str(e)}")
                logger.error(f"Report generation error: {e}", exc_info=True)
    
    st.divider()
    
    # Molecular visualization
    st.markdown("#### 🧬 3D Molecular Structure Viewer")
    
    candidates = st.session_state.scored_candidates[:10]
    
    if candidates:
        drug_names = [c.get('drug_name', 'Unknown') for c in candidates]
        selected_drug = st.selectbox(
            "Select Drug to Visualize",
            drug_names,
            help="View 3D molecular structure of selected candidate"
        )
        
        if selected_drug:
            col1, col2 = st.columns([2, 1])
            
            with col1:
                get_and_display_molecule(selected_drug, width=700, height=500)
            
            with col2:
                selected_candidate = next((c for c in candidates if c.get('drug_name') == selected_drug), None)
                if selected_candidate:
                    st.metric("Score", f"{selected_candidate.get('score', 0)}/100")
                    st.metric("FDA Approval", selected_candidate.get('approval_year', 'Unknown'))
                    st.metric("Rank", f"#{candidates.index(selected_candidate) + 1}")
                    
                    st.markdown("**Pathways:**")
                    for pathway in selected_candidate.get('pathway_targets', []):
                        st.caption(f"• {pathway}")


if __name__ == "__main__":
    main()
