"""
Script to update app.py with reasoning trace functionality
"""

import re

# Read the original app.py
with open('/root/aiCoScientist/app.py', 'r') as f:
    content = f.read()

# 1. Add reasoning_trace to session state initialization
session_state_init = """    defaults = {
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
    }"""

content = re.sub(
    r"defaults = \{[^}]+\}",
    session_state_init,
    content,
    flags=re.DOTALL
)

# 2. Add CSS for reasoning trace
css_addition = """    
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
"""

content = content.replace("</style>", css_addition)

# 3. Add reasoning callback function before run_full_analysis
reasoning_callback_code = '''

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


'''

# Insert before run_full_analysis function
content = content.replace(
    "def run_full_analysis(",
    reasoning_callback_code + "\ndef run_full_analysis("
)

# 4. Update run_full_analysis to use reasoning callback and display live trace
run_analysis_update = '''def run_full_analysis(
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
'''

# Replace the run_full_analysis function
content = re.sub(
    r'def run_full_analysis\([^)]+\):.*?(?=\ndef |\Z)',
    run_analysis_update,
    content,
    flags=re.DOTALL
)

# 5. Add reasoning trace tab to the workbench
# Find the tabs definition and add reasoning trace tab
tabs_pattern = r'tabs = st\.tabs\(\[(.*?)\]\)'
tabs_match = re.search(tabs_pattern, content)

if tabs_match:
    current_tabs = tabs_match.group(1)
    if '"🧠 Reasoning Trace"' not in current_tabs:
        new_tabs = current_tabs.rstrip() + ',\n        "🧠 Reasoning Trace"'
        content = content.replace(
            f'tabs = st.tabs([{current_tabs}])',
            f'tabs = st.tabs([{new_tabs}])'
        )

# Write updated content
with open('/root/aiCoScientist/app.py', 'w') as f:
    f.write(content)

print("✅ Updated app.py with reasoning trace functionality")
print("Added:")
print("  - reasoning_trace to session state")
print("  - CSS styling for terminal-like reasoning display")
print("  - reasoning_callback function for real-time updates")
print("  - render_reasoning_trace_panel function")
print("  - Updated run_full_analysis to use callback")
print("  - Added Reasoning Trace tab to workbench")
