import traceback
from urllib.parse import quote

import pandas as pd
import requests

from protein_metamorphisms_is.sql.model import Accession


def load_accessions_from_csv(self, csv_path, accession_column, csv_tag):
    """
    Loads accessions from a specified CSV file and processes them for data fetching.
    Args:
        csv_path (str): Path to the CSV file.
        accession_column (str): Column name containing accession codes.
        csv_tag (str): Tag associated with the accession for internal tracking.
    """
    try:
        data = pd.read_csv(csv_path)
        accessions = data[accession_column].dropna().unique()
        self.logger.info(f"Loaded {len(accessions)} unique accession codes from CSV.")
        self.process_new_accessions(accessions, csv_tag)
    except Exception as e:
        self.logger.error(f"Failed to load or process CSV: {traceback.format_exc()}")


def fetch_accessions_from_api(self, search_criteria, limit):
    """
    Fetches accession codes from UniProt API based on the specified search criteria.
    Args:
        search_criteria (str): Search parameters for the API query.
        limit (int): Maximum number of accessions to fetch.
    """
    encoded_search_criteria = quote(search_criteria)
    url = f"https://rest.uniprot.org/uniprotkb/stream?query={encoded_search_criteria}&format=list&size={limit}"
    self.logger.info(f"Fetching data from URL: {url}")
    try:
        response = requests.get(url)
        response.raise_for_status()
        accessions = response.text.strip().split("\n")
        self.logger.info(f"Retrieved {len(accessions)} accessions from UniProt API.")
        self.process_new_accessions(accessions, self.conf.get("tag"))
    except requests.RequestException as e:
        self.logger.error(f"Failed to fetch data from UniProt: {e}")


def process_new_accessions(self, accessions, tag):
    """
    Processes newly fetched accession codes, checking against the database to avoid duplicates,
    and saves them if they are new.
    Args:
        accessions (list): List of accession codes to process.
        tag (str): Tag to associate with new accessions.
    """
    self.logger.info(f"Processing {len(accessions)} accessions.")
    existing_accessions = {acc[0] for acc in self.session.query(Accession.accession_code).filter(
        Accession.accession_code.in_(accessions)).all()}
    new_accessions = [Accession(accession_code=acc, primary=True, tag=tag) for acc in accessions if
                      acc not in existing_accessions]
    self.session.bulk_save_objects(new_accessions)
    self.session.commit()
    self.logger.info(f"Added {len(new_accessions)} new accessions to the database.")