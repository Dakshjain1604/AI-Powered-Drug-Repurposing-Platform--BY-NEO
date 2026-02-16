"""
PubMed Literature Retrieval Module
Fetches abstracts from PubMed using Biopython's Entrez API
"""

from Bio import Entrez
import time
import logging
from typing import List, Dict, Optional
from datetime import datetime, timedelta
import json

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class PubMedFetcher:
    """
    Fetch and parse PubMed abstracts using Biopython's Entrez API
    """
    
    def __init__(self, email: str, api_key: Optional[str] = None):
        """
        Initialize PubMed fetcher
        
        Args:
            email: Email address for NCBI (required)
            api_key: NCBI API key (optional, but increases rate limits)
        """
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        
        self.email = email
        self.api_key = api_key
        self.retry_delay = 2  # seconds
        self.max_retries = 3
    
    def fetch_literature(
        self, 
        query: str, 
        max_results: int = 50,
        time_range_years: int = 5,
        start_date: Optional[str] = None,
        end_date: Optional[str] = None
    ) -> List[Dict]:
        """
        Fetch literature abstracts from PubMed
        
        Args:
            query: Search query (e.g., "liver fibrosis TGF-beta")
            max_results: Maximum number of abstracts to retrieve
            time_range_years: Number of years back to search (default: 5)
            start_date: Custom start date (YYYY/MM/DD format)
            end_date: Custom end date (YYYY/MM/DD format)
            
        Returns:
            List of dictionaries containing abstract metadata
        """
        logger.info(f"Fetching literature for query: {query}")
        
        # Calculate date range if not provided
        if not start_date or not end_date:
            end_date_obj = datetime.now()
            start_date_obj = end_date_obj - timedelta(days=365 * time_range_years)
            start_date = start_date_obj.strftime("%Y/%m/%d")
            end_date = end_date_obj.strftime("%Y/%m/%d")
        
        # Build PubMed query with date filter
        full_query = f"{query} AND {start_date}:{end_date}[pdat]"
        logger.info(f"Full query with date filter: {full_query}")
        
        # Search PubMed for PMIDs
        pmids = self._search_pubmed(full_query, max_results)
        
        if not pmids:
            logger.warning("No PMIDs found for query")
            return []
        
        logger.info(f"Found {len(pmids)} PMIDs, fetching details...")
        
        # Fetch detailed records
        abstracts = self._fetch_details(pmids)
        
        logger.info(f"Successfully retrieved {len(abstracts)} abstracts")
        return abstracts
    
    def _search_pubmed(self, query: str, max_results: int) -> List[str]:
        """
        Search PubMed and return list of PMIDs
        
        Args:
            query: Search query
            max_results: Maximum results to return
            
        Returns:
            List of PMIDs
        """
        for attempt in range(self.max_retries):
            try:
                handle = Entrez.esearch(
                    db="pubmed",
                    term=query,
                    retmax=max_results,
                    sort="relevance",
                    retmode="xml"
                )
                record = Entrez.read(handle)
                handle.close()
                
                pmids = record.get("IdList", [])
                return pmids
                
            except Exception as e:
                logger.warning(f"Search attempt {attempt + 1} failed: {e}")
                if attempt < self.max_retries - 1:
                    time.sleep(self.retry_delay)
                else:
                    logger.error("Max retries reached for PubMed search")
                    return []
        
        return []
    
    def _fetch_details(self, pmids: List[str]) -> List[Dict]:
        """
        Fetch detailed records for a list of PMIDs
        
        Args:
            pmids: List of PubMed IDs
            
        Returns:
            List of abstract dictionaries
        """
        abstracts = []
        
        # Fetch in batches to avoid overwhelming the API
        batch_size = 100
        for i in range(0, len(pmids), batch_size):
            batch_pmids = pmids[i:i+batch_size]
            
            for attempt in range(self.max_retries):
                try:
                    # Fetch records
                    handle = Entrez.efetch(
                        db="pubmed",
                        id=batch_pmids,
                        rettype="medline",
                        retmode="xml"
                    )
                    records = Entrez.read(handle)
                    handle.close()
                    
                    # Parse records
                    for record in records.get("PubmedArticle", []):
                        abstract = self._parse_record(record)
                        if abstract:
                            abstracts.append(abstract)
                    
                    # Rate limiting (3 requests/second without API key, 10/second with key)
                    if not self.api_key:
                        time.sleep(0.34)  # ~3 requests per second
                    else:
                        time.sleep(0.11)  # ~10 requests per second
                    
                    break  # Success, exit retry loop
                    
                except Exception as e:
                    logger.warning(f"Fetch attempt {attempt + 1} failed: {e}")
                    if attempt < self.max_retries - 1:
                        time.sleep(self.retry_delay)
                    else:
                        logger.error(f"Max retries reached for batch {i//batch_size + 1}")
        
        return abstracts
    
    def _parse_record(self, record: dict) -> Optional[Dict]:
        """
        Parse a PubMed record into structured format
        
        Args:
            record: Raw PubMed record
            
        Returns:
            Dictionary with abstract metadata or None if parsing fails
        """
        try:
            medline_citation = record.get("MedlineCitation", {})
            article = medline_citation.get("Article", {})
            
            # Extract PMID
            pmid = str(medline_citation.get("PMID", ""))
            
            # Extract title
            title = article.get("ArticleTitle", "")
            
            # Extract authors
            author_list = article.get("AuthorList", [])
            authors = []
            for author in author_list:
                if "LastName" in author and "Initials" in author:
                    authors.append(f"{author['LastName']} {author['Initials']}")
                elif "CollectiveName" in author:
                    authors.append(author["CollectiveName"])
            
            # Extract publication date
            pub_date = article.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {})
            year = pub_date.get("Year", "")
            month = pub_date.get("Month", "01")
            day = pub_date.get("Day", "01")
            
            # Convert month name to number if needed
            month_map = {
                "Jan": "01", "Feb": "02", "Mar": "03", "Apr": "04",
                "May": "05", "Jun": "06", "Jul": "07", "Aug": "08",
                "Sep": "09", "Oct": "10", "Nov": "11", "Dec": "12"
            }
            if month in month_map:
                month = month_map[month]
            
            publication_date = f"{year}-{month.zfill(2)}-{day.zfill(2)}" if year else ""
            
            # Extract journal
            journal = article.get("Journal", {}).get("Title", "")
            
            # Extract DOI
            article_ids = record.get("PubmedData", {}).get("ArticleIdList", [])
            doi = ""
            for article_id in article_ids:
                if article_id.attributes.get("IdType") == "doi":
                    doi = str(article_id)
                    break
            
            # Extract abstract text
            abstract_sections = article.get("Abstract", {}).get("AbstractText", [])
            abstract_text = ""
            if abstract_sections:
                if isinstance(abstract_sections, list):
                    abstract_text = " ".join([str(section) for section in abstract_sections])
                else:
                    abstract_text = str(abstract_sections)
            
            # Extract keywords
            keyword_list = medline_citation.get("KeywordList", [])
            keywords = []
            if keyword_list:
                for kw_group in keyword_list:
                    keywords.extend([str(kw) for kw in kw_group])
            
            # Build result dictionary
            abstract_dict = {
                "pmid": pmid,
                "title": title,
                "authors": authors,
                "publication_date": publication_date,
                "journal": journal,
                "doi": doi,
                "abstract_text": abstract_text,
                "keywords": keywords
            }
            
            return abstract_dict
            
        except Exception as e:
            logger.error(f"Error parsing record: {e}")
            return None
    
    def save_abstracts(self, abstracts: List[Dict], filepath: str):
        """
        Save abstracts to JSON file
        
        Args:
            abstracts: List of abstract dictionaries
            filepath: Path to save file
        """
        try:
            with open(filepath, 'w', encoding='utf-8') as f:
                json.dump(abstracts, f, indent=2, ensure_ascii=False)
            logger.info(f"Saved {len(abstracts)} abstracts to {filepath}")
        except Exception as e:
            logger.error(f"Error saving abstracts: {e}")
    
    def load_abstracts(self, filepath: str) -> List[Dict]:
        """
        Load abstracts from JSON file
        
        Args:
            filepath: Path to JSON file
            
        Returns:
            List of abstract dictionaries
        """
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                abstracts = json.load(f)
            logger.info(f"Loaded {len(abstracts)} abstracts from {filepath}")
            return abstracts
        except Exception as e:
            logger.error(f"Error loading abstracts: {e}")
            raise Exception(f"Failed to load abstracts from {filepath}: {str(e)}")


def fetch_literature_simple(
    query: str,
    email: str = "user@example.com",
    max_results: int = 50,
    time_range_years: int = 5
) -> List[Dict]:
    """
    Convenience function for simple literature fetching
    
    Args:
        query: Search query
        email: Email for NCBI
        max_results: Maximum results
        time_range_years: Years to search back
        
    Returns:
        List of abstract dictionaries
    """
    fetcher = PubMedFetcher(email=email)
    return fetcher.fetch_literature(query, max_results, time_range_years)
