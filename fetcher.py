from typing import List, Dict
from Bio import Entrez

Entrez.email = "jananithirumal08@gmail.com"

def fetch_paper_ids(query: str, max_results: int) -> List[str]:
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    return record["IdList"]

def fetch_paper_details(ids: List[str]) -> List[Dict]:
    handle = Entrez.efetch(db="pubmed", id=",".join(ids), rettype="medline", retmode="xml")
    records = Entrez.read(handle)
    return records["PubmedArticle"]
   
                   
