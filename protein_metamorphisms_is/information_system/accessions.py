from abc import ABC
import traceback
from urllib.parse import quote
import pandas as pd
import requests

from protein_metamorphisms_is.base.task import BaseTaskInitializer
from protein_metamorphisms_is.sql.model import Accession

class AccessionManager(BaseTaskInitializer):
    def __init__(self, conf, session_required=True):
        super().__init__(conf, session_required)

    def load_accessions_from_csv(self):
        """
        Loads accessions from a specified CSV file and processes them for data fetching.
        """
        try:
            csv_path = self.conf['load_accesion_csv']
            accession_column = self.conf['load_accesion_column']
            csv_tag = self.conf['tag']

            data = pd.read_csv(csv_path)
            accessions = data[accession_column].dropna().unique()
            self.logger.info(f"Loaded {len(accessions)} unique accession codes from CSV.")
            self._process_new_accessions(accessions[:50], csv_tag)
        except Exception as e:
            self.logger.error(f"Failed to load or process CSV: {traceback.format_exc()}")

    def fetch_accessions_from_api(self):
        """
        Fetches accession codes from UniProt API based on the specified search criteria.
        """
        try:
            search_criteria = self.conf['search_criteria']
            limit = self.conf['limit']
            tag = self.conf.get('tag')

            encoded_search_criteria = quote(search_criteria)
            url = f"https://rest.uniprot.org/uniprotkb/stream?query={encoded_search_criteria}&format=list&size={limit}"
            self.logger.info(f"Fetching data from URL: {url}")

            response = requests.get(url)
            response.raise_for_status()
            accessions = response.text.strip().split("\n")
            self.logger.info(f"Retrieved {len(accessions)} accessions from UniProt API.")
            self._process_new_accessions(accessions[:50], tag)
        except requests.RequestException as e:
            self.logger.error(f"Failed to fetch data from UniProt: {e}")

    def _process_new_accessions(self, accessions, tag):
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

    def enqueue(self):
        pass

    def store_entry(self, record):
        pass

    def process(self,_):
        pass