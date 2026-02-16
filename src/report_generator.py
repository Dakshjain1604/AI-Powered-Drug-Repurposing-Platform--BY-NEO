"""
PDF Report Generator
Generates comprehensive research reports using ReportLab
"""

from reportlab.lib.pagesizes import letter, A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_JUSTIFY
from reportlab.lib import colors
from datetime import datetime
from typing import List, Dict
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class ReportGenerator:
    """
    Generate comprehensive PDF reports for drug repurposing research
    """
    
    def __init__(self):
        """
        Initialize report generator
        """
        self.styles = getSampleStyleSheet()
        self._setup_custom_styles()
    
    def _setup_custom_styles(self):
        """
        Setup custom paragraph styles
        """
        # Title style
        self.styles.add(ParagraphStyle(
            name='CustomTitle',
            parent=self.styles['Heading1'],
            fontSize=24,
            textColor=colors.HexColor('#1f77b4'),
            spaceAfter=30,
            alignment=TA_CENTER
        ))
        
        # Section heading
        self.styles.add(ParagraphStyle(
            name='SectionHeading',
            parent=self.styles['Heading2'],
            fontSize=16,
            textColor=colors.HexColor('#2c5aa0'),
            spaceAfter=12,
            spaceBefore=12
        ))
        
        # Subsection heading
        self.styles.add(ParagraphStyle(
            name='SubsectionHeading',
            parent=self.styles['Heading3'],
            fontSize=14,
            textColor=colors.HexColor('#4a90d9'),
            spaceAfter=10,
            spaceBefore=10
        ))
    
    def generate_report(
        self,
        filepath: str,
        disease_focus: str,
        time_range: str,
        repurposing_goal: str,
        abstracts: List[Dict],
        candidates: List[Dict],
        hypotheses: List[Dict],
        protocols: List[Dict]
    ) -> bool:
        """
        Generate comprehensive PDF report
        
        Args:
            filepath: Output file path
            disease_focus: Disease being studied
            time_range: Time range for literature search
            repurposing_goal: User's therapeutic goal
            abstracts: Retrieved abstracts
            candidates: Scored candidates
            hypotheses: Generated hypotheses
            protocols: Experimental protocols
            
        Returns:
            True if successful, False otherwise
        """
        try:
            logger.info(f"Generating PDF report: {filepath}")
            
            # Create PDF document
            doc = SimpleDocTemplate(filepath, pagesize=letter,
                                  rightMargin=72, leftMargin=72,
                                  topMargin=72, bottomMargin=18)
            
            # Build story (content)
            story = []
            
            # Title page
            story.extend(self._create_title_page(disease_focus, repurposing_goal))
            story.append(PageBreak())
            
            # Executive summary
            story.extend(self._create_executive_summary(
                disease_focus, len(abstracts), len(candidates), candidates[:3]
            ))
            story.append(Spacer(1, 0.2*inch))
            
            # Disease background
            story.extend(self._create_disease_background(disease_focus))
            story.append(Spacer(1, 0.2*inch))
            
            # Literature search methodology
            story.extend(self._create_literature_methodology(time_range, len(abstracts)))
            story.append(Spacer(1, 0.2*inch))
            
            # Candidate identification and scoring
            story.extend(self._create_candidates_section(candidates[:10]))
            story.append(PageBreak())
            
            # Top 3 candidates with hypotheses
            story.extend(self._create_top_candidates_section(candidates[:3], hypotheses))
            story.append(PageBreak())
            
            # Experimental protocols
            story.extend(self._create_protocols_section(protocols))
            story.append(PageBreak())
            
            # References
            story.extend(self._create_references_section(abstracts[:10]))
            
            # Build PDF
            doc.build(story)
            
            logger.info(f"Report generated successfully: {filepath}")
            return True
            
        except Exception as e:
            logger.error(f"Error generating report: {e}")
            return False
    
    def _create_title_page(self, disease_focus: str, repurposing_goal: str) -> List:
        """
        Create title page elements
        """
        elements = []
        
        # Title
        title = Paragraph("BioScript Drug Repurposing Report", self.styles['CustomTitle'])
        elements.append(title)
        elements.append(Spacer(1, 0.5*inch))
        
        # Disease focus
        disease = Paragraph(f"<b>Disease Focus:</b> {disease_focus}", self.styles['Normal'])
        elements.append(disease)
        elements.append(Spacer(1, 0.2*inch))
        
        # Goal
        goal = Paragraph(f"<b>Repurposing Goal:</b> {repurposing_goal}", self.styles['Normal'])
        elements.append(goal)
        elements.append(Spacer(1, 0.2*inch))
        
        # Date
        date = Paragraph(f"<b>Report Generated:</b> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", 
                        self.styles['Normal'])
        elements.append(date)
        elements.append(Spacer(1, 0.5*inch))
        
        # Disclaimer
        disclaimer = Paragraph(
            "<i>This report is generated by an AI-assisted drug repurposing platform. "
            "All recommendations should be validated by qualified researchers before experimental use.</i>",
            self.styles['Normal']
        )
        elements.append(disclaimer)
        
        return elements
    
    def _create_executive_summary(
        self, 
        disease_focus: str, 
        num_abstracts: int, 
        num_candidates: int,
        top_candidates: List[Dict]
    ) -> List:
        """
        Create executive summary section
        """
        elements = []
        
        elements.append(Paragraph("Executive Summary", self.styles['SectionHeading']))
        
        summary_text = f"""
        This report presents the results of a systematic drug repurposing analysis for {disease_focus}. 
        We analyzed {num_abstracts} recent scientific abstracts from PubMed and identified {num_candidates} 
        FDA-approved drug candidates with potential therapeutic effects through modulation of fibrotic pathways.
        """
        
        elements.append(Paragraph(summary_text, self.styles['Normal']))
        elements.append(Spacer(1, 0.1*inch))
        
        # Top 3 candidates summary
        if top_candidates:
            elements.append(Paragraph("<b>Top 3 Candidates:</b>", self.styles['Normal']))
            elements.append(Spacer(1, 0.1*inch))
            
            for i, candidate in enumerate(top_candidates[:3], 1):
                drug_name = candidate.get('drug_name', 'Unknown')
                score = candidate.get('score', 0)
                pathways = ', '.join(candidate.get('pathway_targets', [])[:2])
                
                candidate_text = f"{i}. <b>{drug_name}</b> (Score: {score}/100) - Targets: {pathways}"
                elements.append(Paragraph(candidate_text, self.styles['Normal']))
                elements.append(Spacer(1, 0.05*inch))
        
        return elements
    
    def _create_disease_background(self, disease_focus: str) -> List:
        """
        Create disease background section
        """
        elements = []
        
        elements.append(Paragraph("Disease Background", self.styles['SectionHeading']))
        
        background_text = f"""
        {disease_focus} is a progressive pathological condition characterized by excessive accumulation 
        of extracellular matrix (ECM) proteins, leading to tissue scarring and organ dysfunction. The primary 
        drivers include chronic inflammation, myofibroblast activation, and dysregulated wound healing responses. 
        Key molecular pathways include TGF-β signaling, SMAD-mediated transcription, and inflammatory cytokine 
        cascades. Drug repurposing offers a promising approach to identify existing FDA-approved medications 
        that can modulate these pathways and potentially slow or reverse fibrotic progression.
        """
        
        elements.append(Paragraph(background_text, self.styles['Normal']))
        
        return elements
    
    def _create_literature_methodology(self, time_range: str, num_abstracts: int) -> List:
        """
        Create literature search methodology section
        """
        elements = []
        
        elements.append(Paragraph("Literature Search Methodology", self.styles['SectionHeading']))
        
        method_text = f"""
        Literature retrieval was performed using the NCBI PubMed database via the Entrez API. 
        The search covered publications from the last {time_range}. A total of {num_abstracts} abstracts 
        were retrieved based on relevance to fibrotic pathways, drug mechanisms, and therapeutic outcomes. 
        Abstract metadata including PMID, title, authors, publication date, and abstract text were extracted 
        and analyzed using natural language processing and LLM-based reasoning.
        """
        
        elements.append(Paragraph(method_text, self.styles['Normal']))
        
        return elements
    
    def _create_candidates_section(self, candidates: List[Dict]) -> List:
        """
        Create candidate identification and scoring section
        """
        elements = []
        
        elements.append(Paragraph("Candidate Identification and Scoring", self.styles['SectionHeading']))
        
        intro_text = """
        Drug candidates were extracted from abstracts and filtered for FDA approval status. 
        Each candidate was scored (0-100) based on: (1) Pathway relevance (40 points), 
        (2) Mechanism alignment with fibrosis biology (30 points), (3) Clinical evidence (20 points), 
        and (4) Market availability (10 points).
        """
        
        elements.append(Paragraph(intro_text, self.styles['Normal']))
        elements.append(Spacer(1, 0.2*inch))
        
        # Create candidates table
        if candidates:
            table_data = [['Rank', 'Drug Name', 'Score', 'Pathways', 'Evidence']]
            
            for i, candidate in enumerate(candidates, 1):
                drug_name = candidate.get('drug_name', 'Unknown')
                score = candidate.get('score', 0)
                pathways = ', '.join(candidate.get('pathway_targets', [])[:2])
                evidence = candidate.get('evidence_type', 'N/A')
                
                table_data.append([str(i), drug_name, str(score), pathways[:30], evidence[:20]])
            
            table = Table(table_data, colWidths=[0.5*inch, 1.5*inch, 0.7*inch, 2*inch, 1.5*inch])
            table.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, 0), 10),
                ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
                ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
                ('GRID', (0, 0), (-1, -1), 1, colors.black)
            ]))
            
            elements.append(table)
        
        return elements
    
    def _create_top_candidates_section(self, top_candidates: List[Dict], hypotheses: List[Dict]) -> List:
        """
        Create top 3 candidates with hypotheses section
        """
        elements = []
        
        elements.append(Paragraph("Top 3 Candidates with Structured Hypotheses", self.styles['SectionHeading']))
        
        for i, (candidate, hypothesis) in enumerate(zip(top_candidates, hypotheses), 1):
            elements.append(Paragraph(f"Candidate #{i}: {candidate.get('drug_name', 'Unknown')}", 
                                    self.styles['SubsectionHeading']))
            
            # Score breakdown
            score_breakdown = candidate.get('score_breakdown', {})
            breakdown_text = f"""
            <b>Total Score:</b> {candidate.get('score', 0)}/100<br/>
            • Pathway Hit: {score_breakdown.get('pathway_hit', 0)}/40<br/>
            • Mechanism Alignment: {score_breakdown.get('mechanism_alignment', 0)}/30<br/>
            • Clinical Evidence: {score_breakdown.get('clinical_evidence', 0)}/20<br/>
            • Market Availability: {score_breakdown.get('market_availability', 0)}/10
            """
            elements.append(Paragraph(breakdown_text, self.styles['Normal']))
            elements.append(Spacer(1, 0.1*inch))
            
            # Hypothesis
            elements.append(Paragraph("<b>Hypothesis:</b>", self.styles['Normal']))
            hypothesis_text = hypothesis.get('hypothesis', 'No hypothesis generated')
            elements.append(Paragraph(hypothesis_text, self.styles['Normal']))
            elements.append(Spacer(1, 0.1*inch))
            
            # Expected outcomes
            outcomes = hypothesis.get('expected_outcomes', [])
            if outcomes:
                elements.append(Paragraph("<b>Expected Outcomes:</b>", self.styles['Normal']))
                for outcome in outcomes:
                    elements.append(Paragraph(f"• {outcome}", self.styles['Normal']))
            
            elements.append(Spacer(1, 0.2*inch))
        
        return elements
    
    def _create_protocols_section(self, protocols: List[Dict]) -> List:
        """
        Create experimental protocols section
        """
        elements = []
        
        elements.append(Paragraph("Experimental Protocols", self.styles['SectionHeading']))
        
        for i, protocol in enumerate(protocols, 1):
            drug_name = protocol.get('drug_name', 'Unknown')
            elements.append(Paragraph(f"Protocol #{i}: {drug_name}", self.styles['SubsectionHeading']))
            
            protocol_text = f"""
            <b>Model System:</b> {protocol.get('model_system', 'N/A')}<br/>
            <b>Assay Type:</b> {protocol.get('assay_type', 'N/A')}<br/>
            <b>Drug Concentration:</b> {protocol.get('drug_concentration_range', 'N/A')}<br/>
            <b>Dosing Schedule:</b> {protocol.get('dosing_schedule', 'N/A')}<br/>
            <b>Timeline:</b> {protocol.get('assay_timeline', 'N/A')}<br/>
            <b>Positive Control:</b> {protocol.get('controls', {}).get('positive', 'N/A')}<br/>
            <b>Negative Control:</b> {protocol.get('controls', {}).get('negative', 'N/A')}<br/>
            <b>Statistical Analysis:</b> {protocol.get('statistical_analysis', 'N/A')}
            """
            
            elements.append(Paragraph(protocol_text, self.styles['Normal']))
            
            # Success metrics
            metrics = protocol.get('success_metrics', [])
            if metrics:
                elements.append(Spacer(1, 0.1*inch))
                elements.append(Paragraph("<b>Success Metrics:</b>", self.styles['Normal']))
                for metric in metrics:
                    elements.append(Paragraph(f"• {metric}", self.styles['Normal']))
            
            elements.append(Spacer(1, 0.2*inch))
        
        return elements
    
    def _create_references_section(self, abstracts: List[Dict]) -> List:
        """
        Create references section
        """
        elements = []
        
        elements.append(Paragraph("Key References", self.styles['SectionHeading']))
        
        for i, abstract in enumerate(abstracts, 1):
            authors = ', '.join(abstract.get('authors', [])[:3])
            if len(abstract.get('authors', [])) > 3:
                authors += ' et al.'
            
            title = abstract.get('title', 'Unknown title')
            journal = abstract.get('journal', 'Unknown journal')
            year = abstract.get('publication_date', '')[:4]
            pmid = abstract.get('pmid', '')
            
            ref_text = f"{i}. {authors} ({year}). {title}. <i>{journal}</i>. PMID: {pmid}"
            elements.append(Paragraph(ref_text, self.styles['Normal']))
            elements.append(Spacer(1, 0.05*inch))
        
        return elements
