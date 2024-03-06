from concurrent.futures import ThreadPoolExecutor, as_completed
from http.client import HTTPException

import requests
from Bio import ExPASy, SwissProt
from sqlalchemy import func, exists
from urllib.parse import quote

from protein_metamorphisms_is.helpers.parser.parser import extract_float
from protein_metamorphisms_is.information_system.base.extractor import ExtractorBase

from protein_metamorphisms_is.sql.model import Accession, Protein, PDBReference


class UniProtExtractor(ExtractorBase):
    """
    A class for extracting and processing data from UniProt, a comprehensive resource for protein sequence and annotation data.
    UniProt provides a rich collection of protein sequence and functional information, which includes protein names,
    descriptions, taxonomic data, and sequence annotations.

    This class extends BioinfoExtractorBase and provides specific implementations for extracting and processing data from UniProt.
    """

    def __init__(self, conf):
        """
        Initialize the UniProtExtractor class.

        Sets up the configuration and logger. Initializes the database session if required.
        Configurations can include number of threads, database connection details, logging settings, and specific
        parameters related to UniProt data extraction.
        """
        super().__init__(conf, session_required=True)

    def start(self):
        """
        Start the data extraction process for UniProt.

        Initiates the process of fetching protein data from UniProt based on predefined search criteria and limits.
        Implements logic for handling the extraction and processing of data in a structured and efficient manner.
        """
        try:
            self.logger.info("Starting UniProt data extraction")
            search_criteria = self.conf.get("search_criteria")
            limit = self.conf.get("limit")
            self.load_access_codes(search_criteria, limit)
            self.extract_entries()

        except Exception as e:
            self.logger.error(f"Error during extraction process: {e}")

    def load_access_codes(self, search_criteria, limit):
        """
        Load access codes from UniProt based on the given search criteria and limit.

        Fetches accession codes from UniProt using RESTful API calls. Accession codes are unique identifiers for protein records in UniProt.
        The function uses these codes to selectively download detailed protein information in later stages.

        Args:
            search_criteria (str): The search criteria for querying UniProt.
            limit (int): The maximum number of results to fetch from UniProt. (A parameter requested by Uniprot with no
            significant impact.)
        """
        encoded_search_criteria = quote(search_criteria)
        url = (
            f"https://rest.uniprot.org/uniprotkb/stream?"
            f"query={encoded_search_criteria}"
            f"&format=list&size={limit}"
        )
        self.logger.info(f"Requested URL: {url}")
        try:
            response = requests.get(url)
            response.raise_for_status()
            accessions = response.text.strip().split("\n")
            self.logger.info(f"Number of Accessions in UniProt: {len(accessions)}")

            found_accessions = []
            new_accessions = []

            for accession_code in accessions:
                exists_query = exists().where(Accession.accession_code == accession_code)
                accession_exists = self.session.query(exists_query).scalar()

                if accession_exists:
                    self.session.query(Accession) \
                        .filter_by(accession_code=accession_code) \
                        .update({"updated_at": func.now(), "disappeared": False})
                    found_accessions.append(accession_code)
                else:
                    new_accession = Accession(accession_code=accession_code, primary=True)
                    self.session.add(new_accession)
                    new_accessions.append(new_accession)

            not_found_accessions = (
                self.session.query(Accession)
                .filter(~Accession.accession_code.in_(accessions))
                .filter(~Accession.disappeared)
                .update({Accession.disappeared: True}, synchronize_session=False)
            )
            self.logger.info(f"Accessions not found: {not_found_accessions}")

            self.session.commit()

            self.logger.info(f"Existing accessions: {len(found_accessions)}")
            self.logger.info(f"New accessions: {len(new_accessions)}")

        except Exception as e:
            self.session.rollback()
            self.logger.error(f"Error: {e}")

    def extract_entries(self):
        """
        Download and process UniProt entries concurrently using multi-threading.

        Uses ThreadPoolExecutor for concurrent downloads, which significantly speeds up the data extraction process,
        especially beneficial when dealing with large datasets.
        """

        self.logger.info("Starting the download of UniProt entries.")

        accessions = self.session.query(Accession).all()
        self.logger.info(f"Total proteins to download: {len(accessions)}")

        max_workers = self.conf.get("max_workers")
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_uniprot_id = {
                executor.submit(
                    self.download_record, accession.accession_code
                ): accession
                for accession in accessions
            }
            for future in as_completed(future_to_uniprot_id):
                uniprot_id = future_to_uniprot_id[future]
                try:
                    data = future.result()
                    if data:
                        self.store_entry(data)
                except Exception as e:
                    self.logger.error(f"Error processing the entry {uniprot_id}: {e}")

    def download_record(self, accession_code):
        """
        Download detailed protein information from UniProt using ExPASy and SwissProt.

        ExPASy is a Bioinformatics Resource Portal which provides access to scientific databases and software tools,
        while SwissProt is a manually annotated and reviewed protein sequence database part of UniProt.

        Args:
            accession_code (str): The accession code of the protein record to download.
        """

        try:
            handle = ExPASy.get_sprot_raw(accession_code)
            record = SwissProt.read(handle)
            return record

        except Exception as e:
            self.logger.error(f"Error downloading the entry {accession_code}: {e}")
            return None

    def store_entry(self, data):
        """
        Stores the downloaded UniProt data in the database.

        Processes and stores detailed protein information including annotations, cross-references, and sequence data.
        The function is designed to handle both new entries and updates to existing records, ensuring data consistency.

        Args:
            data (SwissProt.Record): The UniProt data record to store.
        """
        try:
            exists_query = exists().where(Protein.entry_name == data.entry_name)
            protein_exists = self.session.query(exists_query).scalar()

            if protein_exists:
                protein = self.session.query(Protein).filter_by(entry_name=data.entry_name).first()
            else:
                protein = Protein(entry_name=data.entry_name)
                self.session.add(protein)

            protein.data_class = data.data_class
            protein.molecule_type = data.molecule_type
            protein.sequence_length = data.sequence_length
            protein.created_date = data.created[0]
            protein.sequence = data.sequence
            protein.sequence_update_date = data.sequence_update[0]
            protein.annotation_update_date = data.annotation_update[0]
            protein.description = data.description
            protein.gene_name = str(data.gene_name)
            protein.organism = data.organism
            protein.organelle = data.organelle
            protein.organism_classification = ",".join(
                data.organism_classification
            )
            protein.taxonomy_id = ",".join(data.taxonomy_id)
            protein.host_organism = ",".join(data.host_organism)
            protein.host_taxonomy_id = ",".join(data.host_taxonomy_id)
            protein.comments = "; ".join(data.comments)
            protein.keywords = data.keywords
            protein.protein_existence = data.protein_existence
            protein.seqinfo = data.seqinfo

            self.session.add(protein)

            for accession_code in data.accessions:
                exists_query = exists().where(Accession.accession_code == accession_code)
                accession_exists = self.session.query(exists_query).scalar()

                if not accession_exists:
                    new_accession = Accession(accession_code=accession_code, primary=False)
                    new_accession.protein_entry_name = protein.entry_name
                    self.session.add(new_accession)

            for reference in data.cross_references:
                if reference[0] == "PDB":
                    pdb_ref_exists = self.session.query(exists().where(PDBReference.pdb_id == reference[1])).scalar()
                    if not pdb_ref_exists:
                        pdb_ref = PDBReference(
                            pdb_id=reference[1],
                            method=reference[2],
                            resolution=extract_float(reference[3]),
                            protein=protein,
                        )
                        self.session.add(pdb_ref)

            self.session.commit()

        except HTTPException as e:
            self.logger.error(f"Error while dumping the entry: {e}")
            self.session.rollback()
        finally:
            self.session.close()
