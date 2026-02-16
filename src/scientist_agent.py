"""
Research Agent - Core reasoning engine for biomedical research
Integrates LLM for extraction, scoring, hypothesis generation, and protocol design
Exposes all intermediate data for transparent research workflows
"""

import json
import logging
from typing import List, Dict, Optional, Tuple, Callable
import re
from datetime import datetime
from openai import OpenAI
from src.pathway_database import FIBROTIC_PATHWAYS, get_pathway_weight, search_pathways_by_molecule
from src.fda_database import FDADatabase

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class ScientistAgent:
    """
    Core reasoning engine for biomedical research
    Exposes intermediate datasets for full research transparency
    """
    
    def __init__(
        self, 
        openai_api_key: str,
        model: str = "gpt-4o-mini",
        temperature: float = 0.3,
        reasoning_callback: Optional[Callable] = None
    ):
        """
        Initialize Research Agent
        
        Args:
            openai_api_key: OpenAI API key
            model: LLM model to use
            temperature: Temperature for LLM (lower = more deterministic)
            reasoning_callback: Optional callback function for real-time reasoning updates
                               Signature: callback(step_type, message, data)
        """
        self.client = OpenAI(api_key=openai_api_key)
        self.model = model
        self.temperature = temperature
        self.fda_db = FDADatabase(use_api=False)  # Use local database for speed
        
        # Reasoning trace callback for real-time UI updates
        self.reasoning_callback = reasoning_callback
        self.reasoning_trace = []  # Full reasoning history
        
        # Intermediate data tracking for transparency
        self.all_extracted_drugs = []  # All drugs before FDA filtering
        self.filtered_out_drugs = []  # Drugs filtered out (non-FDA)
        self.extraction_logs = []  # Detailed extraction logs
        self.scoring_logs = []  # Detailed scoring breakdown
        
        self._log_reasoning("initialization", f"✅ Initialized Research Agent with model: {model}", {
            "model": model,
            "temperature": temperature
        })
        logger.info(f"Initialized Research Agent with model: {model}")
    
    def _log_reasoning(self, step_type: str, message: str, data: Dict = None):
        """
        Log reasoning step for transparency
        
        Args:
            step_type: Type of step (e.g., 'extraction', 'scoring', 'filtering')
            message: Human-readable message
            data: Optional data dictionary with details
        """
        trace_entry = {
            "timestamp": datetime.now().isoformat(),
            "step_type": step_type,
            "message": message,
            "data": data or {}
        }
        self.reasoning_trace.append(trace_entry)
        
        # Call callback if provided (for real-time UI updates)
        if self.reasoning_callback:
            try:
                self.reasoning_callback(step_type, message, data)
            except Exception as e:
                logger.warning(f"Callback error: {e}")
    
    def extract_candidates_from_text(
        self, 
        abstracts: List[Dict],
        disease_focus: str = "liver fibrosis"
    ) -> Dict:
        """
        Extract drug candidates from abstracts using LLM
        Returns FULL intermediate data for research transparency
        
        Args:
            abstracts: List of abstract dictionaries
            disease_focus: Disease being studied
            
        Returns:
            Dictionary containing:
              - fda_approved: List of FDA-approved candidates
              - all_extracted: List of all extracted drugs (before filtering)
              - filtered_out: List of non-FDA drugs
              - extraction_logs: Detailed logs
        """
        self._log_reasoning("extraction_start", f"🔍 Starting candidate extraction from {len(abstracts)} abstracts", {
            "total_abstracts": len(abstracts),
            "disease_focus": disease_focus
        })
        logger.info(f"Extracting candidates from {len(abstracts)} abstracts")
        
        # Reset intermediate tracking
        self.all_extracted_drugs = []
        self.filtered_out_drugs = []
        self.extraction_logs = []
        
        # Prepare abstracts text
        self._log_reasoning("llm_preparation", f"📝 Preparing abstracts for LLM analysis (using top 20 for efficiency)", {
            "abstracts_to_process": min(20, len(abstracts))
        })
        abstracts_text = self._prepare_abstracts_for_llm(abstracts[:20])  # Limit for token efficiency
        
        # LLM prompt for extraction
        self._log_reasoning("llm_query", f"🤖 Querying LLM for drug candidate extraction", {
            "model": self.model,
            "temperature": self.temperature,
            "task": "Extract FDA-approved drug candidates with pathway information"
        })
        
        prompt = f"""You are a biomedical research agent analyzing scientific literature for drug repurposing opportunities in {disease_focus}.

Analyze the following abstracts and extract drug candidates mentioned in the context of {disease_focus} or related pathways (TGF-β signaling, fibrosis, ECM remodeling, inflammation).

ABSTRACTS:
{abstracts_text}

For each drug mentioned, extract:
1. Drug name (generic name preferred)
2. Mechanism of action related to fibrosis
3. Pathway targets (e.g., TGF-β signaling, myofibroblast activation, ECM remodeling, inflammation)
4. Evidence type (clinical trial, in-vitro, in-vivo, review)
5. Source PMID

Return ONLY a JSON array of drug objects. Each object should have:
{{
  "drug_name": "string",
  "mechanism": "string",
  "pathway_targets": ["pathway1", "pathway2"],
  "evidence_type": "string",
  "source_pmid": "string",
  "extraction_confidence": 0.0-1.0
}}

Focus on drugs that directly or indirectly affect fibrotic pathways. Return at least 5-10 candidates if present.
"""
        
        try:
            self._log_reasoning("llm_call", "⏳ Waiting for LLM response...", {})
            response = self.client.chat.completions.create(
                model=self.model,
                messages=[
                    {"role": "system", "content": "You are a biomedical research agent specializing in drug repurposing analysis and data extraction."},
                    {"role": "user", "content": prompt}
                ],
                temperature=self.temperature,
                max_tokens=2000
            )
            
            content = response.choices[0].message.content.strip()
            self._log_reasoning("llm_response", "✅ LLM response received", {
                "response_length": len(content),
                "tokens_used": response.usage.total_tokens if hasattr(response, 'usage') else None
            })
            
            # Extract JSON from response
            self._log_reasoning("json_parsing", "📊 Parsing JSON response from LLM", {})
            candidates = self._extract_json_from_response(content)
            
            # Store ALL extracted drugs (before filtering)
            self.all_extracted_drugs = candidates.copy() if candidates else []
            self._log_reasoning("extraction_complete", f"📦 Extracted {len(candidates)} total drug candidates", {
                "total_extracted": len(candidates),
                "candidate_names": [c.get('drug_name', 'Unknown') for c in candidates[:10]]
            })
            
            self.extraction_logs.append({
                "timestamp": datetime.now().isoformat(),
                "total_extracted": len(candidates),
                "source_abstracts": len(abstracts)
            })
            
            # Filter for FDA-approved drugs
            self._log_reasoning("fda_filtering_start", "🔍 Checking FDA approval status for each candidate", {
                "total_to_check": len(candidates)
            })
            
            fda_approved_candidates = []
            for idx, candidate in enumerate(candidates, 1):
                drug_name = candidate.get("drug_name", "").strip()
                
                self._log_reasoning("fda_check", f"Checking FDA status for: {drug_name} ({idx}/{len(candidates)})", {
                    "drug_name": drug_name,
                    "progress": f"{idx}/{len(candidates)}"
                })
                
                # Check FDA approval status
                is_approved = self.fda_db.is_fda_approved(drug_name)
                
                if is_approved:
                    drug_info = self.fda_db.get_drug_info(drug_name)
                    candidate["fda_approved"] = True
                    candidate["approval_year"] = drug_info.get("approval_year") if drug_info else None
                    candidate["indication"] = drug_info.get("indication") if drug_info else "Unknown"
                    fda_approved_candidates.append(candidate)
                    
                    self._log_reasoning("fda_approved", f"✅ FDA-approved: {drug_name}", {
                        "drug_name": drug_name,
                        "approval_year": candidate["approval_year"],
                        "indication": candidate["indication"]
                    })
                    logger.info(f"FDA-approved candidate found: {drug_name}")
                else:
                    # Track filtered-out drugs
                    candidate["fda_approved"] = False
                    candidate["filter_reason"] = "Not FDA approved"
                    self.filtered_out_drugs.append(candidate)
                    
                    self._log_reasoning("fda_rejected", f"❌ Not FDA-approved: {drug_name}", {
                        "drug_name": drug_name,
                        "reason": "Not found in FDA database"
                    })
                    logger.info(f"Non-FDA drug filtered: {drug_name}")
            
            # Update logs
            self.extraction_logs.append({
                "timestamp": datetime.now().isoformat(),
                "fda_approved_count": len(fda_approved_candidates),
                "filtered_out_count": len(self.filtered_out_drugs)
            })
            
            self._log_reasoning("filtering_complete", f"✅ FDA filtering complete: {len(fda_approved_candidates)} approved, {len(self.filtered_out_drugs)} filtered out", {
                "fda_approved_count": len(fda_approved_candidates),
                "filtered_out_count": len(self.filtered_out_drugs),
                "approved_drugs": [c.get('drug_name', 'Unknown') for c in fda_approved_candidates]
            })
            
            logger.info(f"Extracted {len(fda_approved_candidates)} FDA-approved candidates, filtered {len(self.filtered_out_drugs)} non-FDA drugs")
            
            # Return comprehensive results
            return {
                "fda_approved": fda_approved_candidates,
                "all_extracted": self.all_extracted_drugs,
                "filtered_out": self.filtered_out_drugs,
                "extraction_logs": self.extraction_logs
            }
            
        except Exception as e:
            self._log_reasoning("error", f"❌ Error in candidate extraction: {str(e)}", {
                "error": str(e),
                "error_type": type(e).__name__
            })
            logger.error(f"Error in candidate extraction: {e}", exc_info=True)
            raise Exception(f"Failed to extract candidates: {str(e)}. Please check your OpenAI API key and try again.")
    
    def score_candidates(
        self, 
        candidates: List[Dict],
        abstracts: List[Dict]
    ) -> Dict:
        """
        Score candidates based on pathway relevance
        Returns FULL scoring logs for transparency
        Scoring: +40 pathway hit, +30 mechanism alignment, +20 clinical evidence, +10 availability
        
        Args:
            candidates: List of candidate dictionaries
            abstracts: Original abstracts for context
            
        Returns:
            Dictionary containing:
              - scored_candidates: List of candidates with scores
              - scoring_logs: Detailed scoring breakdown for each candidate
        """
        self._log_reasoning("scoring_start", f"📊 Starting candidate scoring for {len(candidates)} candidates", {
            "total_candidates": len(candidates)
        })
        logger.info(f"Scoring {len(candidates)} candidates")
        
        # Reset scoring logs
        self.scoring_logs = []
        scored_candidates = []
        
        for idx, candidate in enumerate(candidates, 1):
            drug_name = candidate.get("drug_name", "Unknown")
            self._log_reasoning("scoring_candidate", f"Scoring {drug_name} ({idx}/{len(candidates)})", {
                "drug_name": drug_name,
                "progress": f"{idx}/{len(candidates)}"
            })
            
            score = 0
            score_breakdown = {}
            
            # 1. Pathway hit score (0-40 points)
            self._log_reasoning("scoring_pathway", f"Calculating pathway score for {drug_name}", {
                "pathway_targets": candidate.get("pathway_targets", [])
            })
            pathway_score = self._calculate_pathway_score(candidate)
            score += pathway_score
            score_breakdown["pathway_hit"] = pathway_score
            
            # 2. Mechanism alignment (0-30 points)
            self._log_reasoning("scoring_mechanism", f"Calculating mechanism alignment for {drug_name}", {
                "mechanism": candidate.get("mechanism", "")[:100]
            })
            mechanism_score = self._calculate_mechanism_score(candidate, abstracts)
            score += mechanism_score
            score_breakdown["mechanism_alignment"] = mechanism_score
            
            # 3. Clinical evidence (0-20 points)
            self._log_reasoning("scoring_clinical", f"Calculating clinical evidence score for {drug_name}", {
                "evidence_type": candidate.get("evidence_type", "")
            })
            clinical_score = self._calculate_clinical_score(candidate)
            score += clinical_score
            score_breakdown["clinical_evidence"] = clinical_score
            
            # 4. Market availability (0-10 points)
            self._log_reasoning("scoring_availability", f"Calculating market availability for {drug_name}", {
                "fda_approved": candidate.get("fda_approved", False),
                "approval_year": candidate.get("approval_year")
            })
            availability_score = self._calculate_availability_score(candidate)
            score += availability_score
            score_breakdown["market_availability"] = availability_score
            
            # Add score to candidate
            candidate["score"] = score
            candidate["score_breakdown"] = score_breakdown
            
            self._log_reasoning("score_complete", f"✅ {drug_name}: Total score = {score}/100", {
                "drug_name": drug_name,
                "total_score": score,
                "breakdown": score_breakdown
            })
            
            # Log detailed scoring
            scoring_log = {
                "drug_name": drug_name,
                "total_score": score,
                "breakdown": score_breakdown,
                "pathway_targets": candidate.get("pathway_targets", []),
                "mechanism": candidate.get("mechanism", ""),
                "evidence_type": candidate.get("evidence_type", "")
            }
            self.scoring_logs.append(scoring_log)
            
            scored_candidates.append(candidate)
            logger.info(f"{candidate['drug_name']}: Total score = {score}")
        
        # Sort by score (descending)
        scored_candidates.sort(key=lambda x: x["score"], reverse=True)
        
        self._log_reasoning("scoring_complete", f"✅ Scoring complete. Top candidate: {scored_candidates[0]['drug_name']} ({scored_candidates[0]['score']} points)", {
            "total_scored": len(scored_candidates),
            "top_3": [{"drug": c["drug_name"], "score": c["score"]} for c in scored_candidates[:3]]
        })
        
        # Return comprehensive results
        return {
            "scored_candidates": scored_candidates,
            "scoring_logs": self.scoring_logs
        }
    
    def _calculate_pathway_score(self, candidate: Dict) -> int:
        """
        Calculate pathway hit score (0-40 points)
        """
        pathway_targets = candidate.get("pathway_targets", [])
        
        if not pathway_targets:
            return 0
        
        # Check for high-value pathways
        high_value_pathways = ["TGF-β signaling", "myofibroblast activation", "SMAD signaling"]
        medium_value_pathways = ["ECM remodeling", "inflammation"]
        
        score = 0
        for pathway in pathway_targets:
            pathway_lower = pathway.lower()
            if any(hvp.lower() in pathway_lower for hvp in high_value_pathways):
                score += 20  # High-value pathway
            elif any(mvp.lower() in pathway_lower for mvp in medium_value_pathways):
                score += 10  # Medium-value pathway
            else:
                score += 5  # Other relevant pathway
        
        return min(score, 40)  # Cap at 40
    
    def _calculate_mechanism_score(self, candidate: Dict, abstracts: List[Dict]) -> int:
        """
        Calculate mechanism alignment score (0-30 points)
        """
        mechanism = candidate.get("mechanism", "").lower()
        
        if not mechanism:
            return 0
        
        # Keywords for fibrosis-relevant mechanisms
        high_relevance = ["fibrosis", "fibrotic", "collagen", "tgf", "myofibroblast", "stellate"]
        medium_relevance = ["inflammation", "inflammatory", "oxidative", "ecm", "matrix"]
        
        score = 0
        for keyword in high_relevance:
            if keyword in mechanism:
                score += 10
        
        for keyword in medium_relevance:
            if keyword in mechanism:
                score += 5
        
        return min(score, 30)  # Cap at 30
    
    def _calculate_clinical_score(self, candidate: Dict) -> int:
        """
        Calculate clinical evidence score (0-20 points)
        """
        evidence_type = candidate.get("evidence_type", "").lower()
        
        if "clinical trial" in evidence_type or "clinical study" in evidence_type:
            return 20
        elif "in-vivo" in evidence_type or "animal" in evidence_type:
            return 15
        elif "in-vitro" in evidence_type or "cell" in evidence_type:
            return 10
        elif "review" in evidence_type or "meta-analysis" in evidence_type:
            return 12
        else:
            return 5
    
    def _calculate_availability_score(self, candidate: Dict) -> int:
        """
        Calculate market availability score (0-10 points)
        """
        if candidate.get("fda_approved", False):
            approval_year = candidate.get("approval_year")
            if approval_year and approval_year < 2000:
                return 10  # Generic likely available
            else:
                return 8  # Newer drug, may be more expensive
        
        return 0
    
    def generate_structured_hypothesis(
        self, 
        top_candidates: List[Dict],
        disease_focus: str = "liver fibrosis"
    ) -> List[Dict]:
        """
        Generate structured scientific hypotheses for top candidates
        
        Args:
            top_candidates: Top 3 candidates (sorted by score)
            disease_focus: Disease being studied
            
        Returns:
            List of hypothesis dictionaries
        """
        self._log_reasoning("hypothesis_start", f"💡 Generating hypotheses for top {len(top_candidates[:3])} candidates", {
            "candidates": [c.get("drug_name", "Unknown") for c in top_candidates[:3]]
        })
        logger.info(f"Generating hypotheses for {len(top_candidates)} candidates")
        
        hypotheses = []
        
        for rank, candidate in enumerate(top_candidates[:3], start=1):
            drug_name = candidate.get("drug_name", "Unknown")
            pathways = candidate.get("pathway_targets", [])
            mechanism = candidate.get("mechanism", "Unknown mechanism")
            
            self._log_reasoning("hypothesis_generation", f"Generating hypothesis for {drug_name} (Rank #{rank})", {
                "drug_name": drug_name,
                "rank": rank,
                "pathways": pathways
            })
            
            # Generate hypothesis using template and LLM
            prompt = f"""Generate a structured scientific hypothesis for drug repurposing research.

DRUG: {drug_name}
MECHANISM: {mechanism}
PATHWAYS: {', '.join(pathways)}
DISEASE: {disease_focus}

Generate a hypothesis following this structure:
"Because [drug X] modulates [pathway Y] through [mechanism], it will reduce [specific outcome Z] in {disease_focus}."

Also provide:
1. Expected biological outcomes (2-3 specific outcomes, e.g., "reduced α-SMA expression", "decreased collagen-I deposition")
2. Reasoning chain (2-3 sentences explaining the biological logic)

Return as JSON:
{{
  "hypothesis": "complete hypothesis statement",
  "expected_outcomes": ["outcome1", "outcome2", "outcome3"],
  "reasoning": "biological reasoning explanation",
  "confidence_score": 0.0-1.0
}}
"""
            
            try:
                response = self.client.chat.completions.create(
                    model=self.model,
                    messages=[
                        {"role": "system", "content": "You are a biomedical research agent generating testable, evidence-based hypotheses."},
                        {"role": "user", "content": prompt}
                    ],
                    temperature=self.temperature,
                    max_tokens=500
                )
                
                content = response.choices[0].message.content.strip()
                hypothesis_data = self._extract_json_from_response(content)
                
                if hypothesis_data:
                    hypothesis_dict = {
                        "candidate_rank": rank,
                        "drug_name": drug_name,
                        "hypothesis": hypothesis_data.get("hypothesis", ""),
                        "supporting_literature": [{"pmid": candidate.get("source_pmid", ""), "excerpt": mechanism}],
                        "expected_outcomes": hypothesis_data.get("expected_outcomes", []),
                        "reasoning": hypothesis_data.get("reasoning", ""),
                        "confidence_score": hypothesis_data.get("confidence_score", 0.7)
                    }
                    hypotheses.append(hypothesis_dict)
                    
                    self._log_reasoning("hypothesis_complete", f"✅ Hypothesis generated for {drug_name}", {
                        "drug_name": drug_name,
                        "hypothesis": hypothesis_data.get("hypothesis", "")[:100] + "..."
                    })
                    logger.info(f"Generated hypothesis for {drug_name}")
                
            except Exception as e:
                self._log_reasoning("hypothesis_error", f"❌ Error generating hypothesis for {drug_name}: {str(e)}", {
                    "drug_name": drug_name,
                    "error": str(e)
                })
                logger.error(f"Error generating hypothesis for {drug_name}: {e}")
        
        self._log_reasoning("hypothesis_generation_complete", f"✅ Generated {len(hypotheses)} hypotheses", {
            "total_hypotheses": len(hypotheses)
        })
        
        return hypotheses
    
    def design_experiment_protocol(
        self, 
        candidate: Dict,
        cell_model: str = "LX-2 human hepatic stellate cells"
    ) -> Dict:
        """
        Design detailed in-vitro experimental protocol
        
        Args:
            candidate: Candidate drug dictionary
            cell_model: Cell model to use
            
        Returns:
            Dictionary with experimental protocol
        """
        drug_name = candidate.get("drug_name", "Unknown")
        mechanism = candidate.get("mechanism", "Unknown")
        pathways = candidate.get("pathway_targets", [])
        
        self._log_reasoning("protocol_design", f"🧪 Designing experimental protocol for {drug_name}", {
            "drug_name": drug_name,
            "cell_model": cell_model
        })
        logger.info(f"Designing experiment protocol for {drug_name}")
        
        prompt = f"""Design a detailed in-vitro experimental protocol for testing a drug repurposing candidate.

DRUG: {drug_name}
MECHANISM: {mechanism}
TARGET PATHWAYS: {', '.join(pathways)}
CELL MODEL: {cell_model}
DISEASE: liver fibrosis

Design a protocol including:
1. Primary assay type (Western blot, ELISA, qPCR, or immunofluorescence)
2. Drug concentration range (e.g., "1-100 µM")
3. Dosing schedule (e.g., "Single dose at time 0, harvest at 24h, 48h, 72h")
4. Assay timeline (specific timepoints)
5. Success metrics (e.g., "≥30% reduction in α-SMA", "≥30% reduction in Collagen-I")
6. Positive control (e.g., "TGF-β1 (10 ng/mL)")
7. Negative control (e.g., "vehicle control (DMSO)")
8. Statistical analysis plan (e.g., "Student's t-test or ANOVA with Bonferroni correction, n=3 replicates")

Return as JSON:
{{
  "assay_type": "string",
  "drug_concentration_range": "string",
  "dosing_schedule": "string",
  "assay_timeline": "string",
  "success_metrics": ["metric1", "metric2"],
  "positive_control": "string",
  "negative_control": "string",
  "statistical_analysis": "string"
}}
"""
        
        try:
            response = self.client.chat.completions.create(
                model=self.model,
                messages=[
                    {"role": "system", "content": "You are a research protocol designer specializing in cell biology experimental design."},
                    {"role": "user", "content": prompt}
                ],
                temperature=self.temperature,
                max_tokens=600
            )
            
            content = response.choices[0].message.content.strip()
            protocol_data = self._extract_json_from_response(content)
            
            if protocol_data:
                protocol = {
                    "drug_name": drug_name,
                    "model_system": cell_model,
                    "assay_type": protocol_data.get("assay_type", "Western blot"),
                    "drug_concentration_range": protocol_data.get("drug_concentration_range", "1-100 µM"),
                    "dosing_schedule": protocol_data.get("dosing_schedule", "Single dose at 0h"),
                    "assay_timeline": protocol_data.get("assay_timeline", "0h, 24h, 48h, 72h"),
                    "success_metrics": protocol_data.get("success_metrics", ["≥30% reduction in α-SMA"]),
                    "controls": {
                        "positive": protocol_data.get("positive_control", "TGF-β1 (10 ng/mL)"),
                        "negative": protocol_data.get("negative_control", "vehicle control (DMSO)")
                    },
                    "statistical_analysis": protocol_data.get("statistical_analysis", "Student's t-test, n=3")
                }
                
                self._log_reasoning("protocol_complete", f"✅ Protocol designed for {drug_name}", {
                    "drug_name": drug_name,
                    "assay_type": protocol["assay_type"]
                })
                logger.info(f"Protocol designed for {drug_name}")
                return protocol
            
        except Exception as e:
            self._log_reasoning("protocol_error", f"❌ Error designing protocol for {drug_name}: {str(e)}", {
                "drug_name": drug_name,
                "error": str(e)
            })
            logger.error(f"Error designing protocol for {drug_name}: {e}")
        
        # Return default protocol if LLM fails
        return self._get_default_protocol(drug_name, cell_model)
    
    def _prepare_abstracts_for_llm(self, abstracts: List[Dict]) -> str:
        """
        Prepare abstracts text for LLM input
        """
        text_parts = []
        for i, abstract in enumerate(abstracts, 1):
            pmid = abstract.get("pmid", "Unknown")
            title = abstract.get("title", "")
            abstract_text = abstract.get("abstract_text", "")
            
            text_parts.append(f"[{i}] PMID: {pmid}\nTitle: {title}\nAbstract: {abstract_text[:500]}...\n")
        
        return "\n".join(text_parts)
    
    def _extract_json_from_response(self, response: str) -> Optional[Dict]:
        """
        Extract JSON from LLM response (handles markdown code blocks)
        """
        try:
            # Try direct JSON parse
            return json.loads(response)
        except:
            # Try to extract from code block
            json_match = re.search(r'```(?:json)?\s*(\[.*?\]|\{.*?\})\s*```', response, re.DOTALL)
            if json_match:
                try:
                    return json.loads(json_match.group(1))
                except:
                    pass
            
            # Try to find JSON array or object
            json_match = re.search(r'(\[.*?\]|\{.*?\})', response, re.DOTALL)
            if json_match:
                try:
                    return json.loads(json_match.group(1))
                except:
                    pass
        
        logger.warning("Could not extract JSON from LLM response")
        return None
    
    def _get_default_protocol(self, drug_name: str, cell_model: str) -> Dict:
        """
        Return default experimental protocol template
        """
        return {
            "drug_name": drug_name,
            "model_system": cell_model,
            "assay_type": "Western blot",
            "drug_concentration_range": "1-100 µM",
            "dosing_schedule": "Single dose at time 0h, harvest at indicated timepoints",
            "assay_timeline": "0h, 24h, 48h, 72h",
            "success_metrics": ["≥30% reduction in α-SMA", "≥30% reduction in Collagen-I"],
            "controls": {
                "positive": "TGF-β1 (10 ng/mL) to induce fibrosis",
                "negative": "vehicle control (DMSO 0.1%)"
            },
            "statistical_analysis": "Student's t-test or one-way ANOVA with Bonferroni correction, n=3 biological replicates, p<0.05"
        }
