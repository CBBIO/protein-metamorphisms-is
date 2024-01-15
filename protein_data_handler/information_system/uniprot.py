from asyncio import as_completed
from concurrent.futures import ThreadPoolExecutor, as_completed
from http.client import HTTPException

import requests
from Bio import ExPASy, SwissProt
from sqlalchemy import func, exists
from sqlalchemy.exc import NoResultFound
from urllib.parse import quote

from protein_data_handler.helpers.config.yaml import read_yaml_config
from protein_data_handler.helpers.parser.parser import extract_float, process_chain_string
from protein_data_handler.information_system.base.bioinfo_extractor import BioinfoExtractorBase

from protein_data_handler.sql.model import Accession, Protein, PDBReference, UniprotChains, GOTerm


class UniProtExtractor(BioinfoExtractorBase):
    """
    A class for extracting data from UniProt.

    This class extends BioinfoExtractorBase and provides specific implementations
    for extracting and processing data from UniProt.

    Args:
        conf (dict): Configuration dictionary containing necessary parameters.
    """

    def __init__(self, conf):
        """
        Initialize the UniProtExtractor class.

        Sets up the configuration and logger. Initializes the database session if required.
        """
        super().__init__(conf, session_required=True)

    def start(self):
        """
        Start the data extraction process for UniProt.

        Overrides the abstract method in BioinfoExtractorBase and
        implements the specific logic for extracting data from UniProt.
        """
        try:
            self.logger.info("Starting UniProt data extraction")
            # Define your search criteria and limit
            search_criteria = self.conf.get("search_criteria")
            limit = self.conf.get("limit")
            # Load access codes and extract entries
            self.load_access_codes(search_criteria, limit)
            self.extract_entries(self.conf['max_workers'])

        except Exception as e:
            self.logger.error(f"Error during extraction process: {e}")

    def load_access_codes(self, search_criteria, limit):
        """
        Load access codes from UniProt based on the given search criteria and limit.

        Args:
            search_criteria (str): The search criteria for querying UniProt.
            limit (int): The maximum number of results to fetch from UniProt.
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
                    # Actualizar el registro existente
                    self.session.query(Accession) \
                        .filter_by(accession_code=accession_code) \
                        .update({"updated_at": func.now(), "disappeared": False})
                    found_accessions.append(accession_code)
                else:
                    # Crear un nuevo registro
                    new_accession = Accession(accession_code=accession_code, primary=True)
                    self.session.add(new_accession)
                    new_accessions.append(new_accession)

            # Actualizar los códigos de acceso que no se encontraron en la búsqueda actual
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

    def extract_entries(self, max_workers=10):
        """
        Download and process UniProt entries concurrently using multi-threading.

        Retrieves access codes from the database and downloads detailed information
        from UniProt concurrently. Uses `ThreadPoolExecutor` to handle multiple simultaneous downloads.

        :param session: SQLAlchemy session for the database.
        :type session: sqlalchemy.orm.session.Session
        :param max_workers: Maximum number of threads for downloads.
        :type max_workers: int, optional
        :raises Exception: For failures in download or processing operations.
        """

        self.logger.info("Starting the download of UniProt entries.")

        accessions = self.session.query(Accession).all()
        self.logger.info(f"Total proteins to download: {len(accessions)}")

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
                    # Get the exception traceback information
                    self.logger.error(f"Error processing the entry {uniprot_id}: {e}")

    def download_record(self, accession_code):
        """
        Download information about a protein from ExPASy/SwissProt.

        :param accession_code: Access code for the download.
        :type accession_code: str

        :raises Exception: If there are errors in the download or processing.

        :return: The protein record or None in case of error.
        :rtype: SwissProt.Record or None
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
        Stores UniProt data in the database.

        Takes a UniProt record and stores it in the database.
        Updates existing entries with the same name or creates a new one
        if it does not exist. Also manages cross-references for each
        UniProt entry.

        :param data: Data of the UniProt entry.
        :type data: SwissProt.Record
        :param session: SQLAlchemy session for database operations.
        :type session: sqlalchemy.orm.session.Session
        :raises Exception: For errors in database operations.
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
                    chains = reference[4].split(',')
                    for chain_obj in chains:
                        chain_name, start, end = process_chain_string(chain_obj)
                        pdb_id = (self.session.query(PDBReference)
                                  .filter_by(pdb_id=reference[1]).first().id)

                        chain = (
                            self.session.query(UniprotChains)
                            .filter_by(pdb_reference_id=pdb_id, chain=chain_name)
                            .first()
                        )
                        if chain is None:
                            chain = UniprotChains(pdb_reference_id=pdb_id,
                                                  chain=chain_name,
                                                  seq_start=start,
                                                  seq_end=end)
                        chain.insert_sequence(protein.sequence)
                        self.session.add(chain)
                elif reference[0] == "GO":
                    go_term_exists = self.session.query(exists().where(GOTerm.go_id == reference[1])).scalar()
                    if not go_term_exists:
                        go_term = GOTerm(
                            go_id=reference[1],
                            category=reference[2].split(":")[0],
                            description=reference[2].split(":")[1],
                            protein=protein,
                        )
                        self.session.add(go_term)

            self.session.commit()

        except HTTPException as e:
            self.logger.error(f"Error while dumping the entry: {e}")
            self.session.rollback()
        finally:
            self.session.close()



