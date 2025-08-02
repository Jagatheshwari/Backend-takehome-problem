Python Program (get-papers-list.py)
import argparse
import csv
from Bio import Entrez
import sys
import xml.etree.ElementTree as ET

# Set your email to comply with NCBI's usage policy
Entrez.email = "jananithirumal08@gmail.com"

# A basic list of keywords to identify non-academic affiliations.
# This can be expanded for more comprehensive filtering.
PHARMA_KEYWORDS = {
   "pharmaceutical", "biotech", "inc.", "gmbh", "ltd.", "llc", "corp.",
   "company", "laboratories", "ag", "s.a.", "plc", "s.p.a",
   # Add specific company names for better accuracy
   "pfizer", "novartis", "roche", "merck", "gsk", "bayer",
   "johnson & johnson", "amgen", "genentech", "astrazeneca",
}

def is_pharma_affiliation(affiliation):
   """Checks if an affiliation string contains any pharma/biotech keywords."""
   if not affiliation:
       return False
   lower_affiliation = affiliation.lower()
   return any(keyword in lower_affiliation for keyword in PHARMA_KEYWORDS)

def fetch_papers(query, debug=False):
   """
   Fetches papers from PubMed based on the query and filters them
   for non-academic affiliations.
   """
   try:
       if debug:
           print(f"Searching PubMed for query: '{query}'...")

       handle = Entrez.esearch(db="pubmed", term=query, retmax=1000)
       record = Entrez.read(handle)
       handle.close()

       id_list = record['IdList']
       if not id_list:
           print("No papers found for the given query.")
           return []

       if debug:
           print(f"Found {len(id_list)} papers. Fetching details...")

       handle = Entrez.efetch(db="pubmed", id=id_list, rettype="fasta", retmode="xml")
       xml_data = handle.read()
       handle.close()

       root = ET.fromstring(xml_data)
       papers = []

       for article in root.findall(".//PubmedArticle"):
           pubmed_id = article.findtext(".//PMID")
           title = article.findtext(".//ArticleTitle")
           pub_date_element = article.find(".//PubDate")
           publication_date = ""
           if pub_date_element is not None:
               year = pub_date_element.findtext("Year")
               month = pub_date_element.findtext("Month")
               day = pub_date_element.findtext("Day")
               publication_date = f"{year}-{month}-{day}" if year else ""
           
           non_academic_authors = []
           company_affiliations = []
           corresponding_author_email = ""
           
           author_list = article.find(".//AuthorList")
           if author_list is not None:
               for author in author_list.findall(".//Author"):
                   affiliation_info = author.findtext(".//Affiliation")
                   if affiliation_info and is_pharma_affiliation(affiliation_info):
                       # Construct author name
                       last_name = author.findtext(".//LastName")
                       fore_name = author.findtext(".//ForeName")
                       author_name = f"{fore_name} {last_name}".strip()
                       if author_name:
                           non_academic_authors.append(author_name)
                       
                       # Extract company names (simple approach)
                       for keyword in PHARMA_KEYWORDS:
                           if keyword in affiliation_info.lower() and keyword not in ["inc.", "ltd."]: # Avoid generic keywords
                               company_affiliations.append(affiliation_info)
                               break
           
           # Extract corresponding author email (PubMed XML schema can be complex here)
           # This is a best-effort attempt as the email is not always present or in a standard location.
           for author in author_list.findall(".//Author"):
               for affiliation_info in author.findall(".//AffiliationInfo/Affiliation"):
                   if "Corresponding Author" in affiliation_info.text:
                       # Email is often in a subsequent block
                       for identifier in author.findall(".//Identifier"):
                           if identifier.get("Source") == "E-mail":
                               corresponding_author_email = identifier.text
                               break
                   if corresponding_author_email:
                       break
               if corresponding_author_email:
                   break
           
           if non_academic_authors:
               papers.append({
                   "PubmedID": pubmed_id,
                   "Title": title,
                   "Publication Date": publication_date,
                   "Non-academic Author(s)": ", ".join(non_academic_authors),
                   "Company Affiliation(s)": ", ".join(set(company_affiliations)), # use set to avoid duplicates
                   "Corresponding Author Email": corresponding_author_email,
               })
       
       return papers
   
   except Exception as e:
       print(f"An error occurred: {e}", file=sys.stderr)
       return []

def main():
   parser = argparse.ArgumentParser(
       description="Fetch research papers from PubMed with non-academic authors."
   )
   parser.add_argument("query", help="The PubMed search query.")
   parser.add_argument("-d", "--debug", action="store_true", help="Print debug information.")
   parser.add_argument("-f", "--file", type=str, help="Specify filename to save results as CSV.")

   args = parser.parse_args()

   results = fetch_papers(args.query, args.debug)

   if not results:
       return

   fieldnames = [
       "PubmedID", "Title", "Publication Date",
       "Non-academic Author(s)", "Company Affiliation(s)",
       "Corresponding Author Email"
   ]

   if args.file:
       try:
           with open(args.file, 'w', newline='', encoding='utf-8') as csvfile:
               writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
               writer.writeheader()
               writer.writerows(results)
           print(f"Successfully saved {len(results)} papers to '{args.file}'.")
       except IOError as e:
           print(f"Error writing to file '{args.file}': {e}", file=sys.stderr)
   else:
       print("--- Results ---")
       for paper in results:
           for key, value in paper.items():
               print(f"{key}: {value}")
           print("-" * 20)

if __name__ == "__main__":
   main()
