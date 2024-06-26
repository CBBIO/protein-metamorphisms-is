import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from http.client import HTTPException

import pandas as pd
import requests
from Bio import ExPASy, SwissProt
from sqlalchemy import func, exists
from urllib.parse import quote

from protein_metamorphisms_is.helpers.logger.logger import setup_logger
from protein_metamorphisms_is.helpers.parser.parser import extract_float, process_chain_string
from protein_metamorphisms_is.information_system.base.extractor import ExtractorBase
from protein_metamorphisms_is.sql.model import Accession, ProteinGOTermAssociation, GOTerm, PDBReference, Sequence, \
    Protein


class UniProtExtractor(ExtractorBase):
    """
    Handles the extraction and processing of data from UniProt. This includes retrieving protein sequences,
    annotations, and functional data. The class extends ExtractorBase to utilize base functionalities and ensures
    integration with a database to store and manage the extracted data.
    """

    def __init__(self, conf):
        """
        Initializes the extractor with configuration for logging and database access.
        Args:
            conf (dict): Configuration dictionary specifying operational parameters.
        """
        super().__init__(conf, session_required=True)
        self.logger.info("UniProtExtractor initialized with configuration.")

    def set_targets(self):
        """
        Determines the data extraction targets from either a CSV file or a live API based on provided configuration.
        """
        csv_path = self.conf.get("load_accesion_csv")
        accession_column = self.conf.get("load_accesion_column")
        search_criteria = self.conf.get("search_criteria")
        limit = self.conf.get("limit")
        csv_tag = self.conf.get("csv_tag")

        if csv_path:
            self._load_access_from_csv(csv_path, accession_column, csv_tag)
        else:
            self.logger.warning("CSV path not provided; attempting to fetch accessions using API.")
            if search_criteria:
                self._fetch_accessions_from_api(search_criteria, limit)
            else:
                self.logger.error("No valid search criteria provided. Data extraction cannot proceed.")

    def _load_access_from_csv(self, csv_path, accession_column, csv_tag):
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
            self._process_new_accessions(accessions, csv_tag)
        except Exception as e:
            self.logger.error(f"Failed to load or process CSV: {traceback.format_exc()}")

    def _fetch_accessions_from_api(self, search_criteria, limit):
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
            self._process_new_accessions(accessions, self.conf.get("tag"))
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

    def fetch(self):
        """
        Initiates the download of UniProt entries. Downloads are processed concurrently using ThreadPoolExecutor
        to manage multiple downloads efficiently. It handles exceptions and logs appropriate errors or warnings.
        """
        self.logger.info("Starting the download of UniProt entries.")

        accessions = self.session.query(Accession).all()
        self.logger.info(f"Total proteins to download: {len(accessions)}")

        max_workers = self.conf.get("max_workers", 10)  # Default to 10 if not specified in config
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_accession = {
                executor.submit(self._download_record, accession.accession_code): accession
                for accession in accessions
            }
            for future in as_completed(future_to_accession):
                accession = future_to_accession[future]
                try:
                    data = future.result()
                    if data:
                        self.data_queue.put(data)  # Add the fetched data to the queue
                        self.logger.info(f"Record for accession {accession.accession_code} added to the queue.")
                    else:
                        self.logger.warning(f"No data found for accession {accession.accession_code}")
                except Exception as e:
                    self.logger.error(f"Error processing the entry for accession {accession.accession_code}: {e}")

        # Signal to the consuming task(s) that fetching is complete by adding None or a similar sentinel value
        self.logger.info("All data fetching and queuing completed.")

    def _download_record(self, accession_code):
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
            self.logger.debug(f"Attempting to store or update protein data for entry_name {data.entry_name}.")
            protein = self._get_or_create_protein(data)
            self._update_protein_details(protein, data)
            self._handle_accessions(protein, data.accessions)
            self._handle_cross_references(protein, data.cross_references)
            self.session.commit()
            self.logger.info(f"Successfully stored or updated data for protein {data.entry_name}.")
        except HTTPException as e:
            self.logger.error(f"HTTP error occurred while processing {data.entry_name}: {e}")
            self.session.rollback()
        except Exception as e:
            self.logger.error(f"An unexpected error occurred while processing {data.entry_name}: {e}")
            self.session.rollback()

    def _get_or_create_protein(self, data):
        """
        Retrieves an existing protein record from the database using the entry_name from data or creates a new one if it does not exist.
        Args:
            data (SwissProt.Record): Data containing the entry_name of the protein.
        Returns:
            Protein: The retrieved or newly created protein object.
        """
        protein = self.session.query(Protein).filter_by(entry_name=data.entry_name).first()
        if not protein:
            protein = Protein(entry_name=data.entry_name)
            self.session.add(protein)
        return protein

    def _update_protein_details(self, protein, data):
        """
        Updates the protein details based on the provided SwissProt data. All related protein attributes like sequence,
        molecule type, organism details, etc., are updated.
        Args:
            protein (Protein): The protein object to update.
            data (SwissProt.Record): Data used to update the protein object.
        """
        protein.sequence = self._get_or_create_sequence(data.sequence)
        protein.data_class = data.data_class
        protein.molecule_type = data.molecule_type
        protein.sequence_length = data.sequence_length
        protein.created_date = data.created[0]
        protein.sequence_update_date = data.sequence_update[0]
        protein.annotation_update_date = data.annotation_update[0]
        protein.description = data.description
        protein.gene_name = str(data.gene_name)
        protein.organism = data.organism
        protein.organelle = data.organelle
        protein.organism_classification = ",".join(data.organism_classification)
        protein.taxonomy_id = ",".join(data.taxonomy_id)
        protein.host_organism = ",".join(data.host_organism)
        protein.host_taxonomy_id = ",".join(data.host_taxonomy_id)
        protein.comments = "; ".join(data.comments)
        protein.keywords = data.keywords
        protein.protein_existence = data.protein_existence
        protein.seqinfo = data.seqinfo
        self.session.add(protein)

    def _get_or_create_sequence(self, sequence):
        """
        Retrieves or creates a sequence entity in the database.
        Args:
            sequence (str): Amino acid sequence of a protein.
        Returns:
            Sequence: The retrieved or newly created sequence object.
        """
        try:
            existing_sequence = self.session.query(Sequence).filter_by(sequence=sequence).first()
            if not existing_sequence:
                existing_sequence = Sequence(sequence=sequence)
                self.session.add(existing_sequence)
            return existing_sequence
        except Exception as e:
            self.logger.error(f"Failed to retrieve or create sequence: {e}")
            raise

    def _handle_accessions(self, protein, accessions):
        """
        Processes and links accession codes to the specified protein, ensuring no duplicates in the database.
        Args:
            protein (Protein): The protein object to which the accession codes will be linked.
            accessions (list): List of accession codes.
        """
        try:
            for accession_code in accessions:
                accession = self._get_or_create_accession(accession_code)
                accession.protein_entry_name = protein.entry_name
                self.session.add(accession)
        except Exception as e:
            self.logger.error(f"Error handling accession {accession_code} from {accessions}: {e}")
            raise

    def _get_or_create_accession(self, accession_code):
        """
        Checks if an association between a protein and a GO term exists in the database and creates it if not.
        Args:
            entry_name (str): The name of the protein entry.
            go_id (str): The GO term identifier.
        Returns:
            ProteinGOTermAssociation: The newly created association, or None if it already exists.
        """
        try:
            exists_query = exists().where(Accession.accession_code == accession_code)
            accession_exists = self.session.query(exists_query).scalar()
            if accession_exists:
                return self.session.query(Accession).filter(Accession.accession_code == accession_code).first()
            else:
                return Accession(accession_code=accession_code, primary=False)
        except Exception as e:
            self.logger.error(f"Failed to retrieve or create accession {accession_code}: {e}")
            raise

    def _handle_cross_references(self, protein, cross_references):
        """
        Manages the cross-references associated with the protein, such as database links to PDB and GO terms.
        Args:
            protein (Protein): The protein object to manage.
            cross_references (list): List of cross-reference data.
        """
        try:
            for reference in cross_references:
                if reference[0] == "PDB":
                    self._handle_pdb_reference(protein, reference)
                elif reference[0] == "GO":
                    self._handle_go_reference(protein, reference)
        except Exception as e:
            self.logger.error(f"Failed to handle cross references for {protein.entry_name}: {e}")
            raise

    def _handle_pdb_reference(self, protein, reference):
        """
        Specific handler for PDB references, extracting relevant segment if specified and storing it.
        Args:
            protein (Protein): The protein object associated with the PDB reference.
            reference (list): PDB reference data, including sequence positions and PDB ID.
        """
        try:
            if reference[0] == "PDB":
                segments = reference[4].split(", ")
                for segment in segments:
                    chain_name, start, end = process_chain_string(segment)
                    if start is not None and end is not None:
                        sequence = protein.sequence.sequence[start - 1:end]  # Adjust indices for Python (0-based index)
                        pdb_ref = self._get_or_create_pdb_reference(reference, sequence)
                        pdb_ref.protein = protein
                        self.session.add(pdb_ref)
        except Exception as e:
            self.logger.error(f"Error handling PDB reference {reference}: {e}")
            raise

    def _get_or_create_pdb_reference(self, reference, sequence):
        """
        Retrieves or creates a PDB reference in the database based on the provided PDB ID.
        If the PDB reference does not exist, it creates a new one using the associated sequence.
        Args:
            reference (list): PDB reference data including PDB ID, method, and resolution.
            sequence (str): Amino acid sequence associated with the PDB entry.
        Returns:
            PDBReference: The retrieved or newly created PDBReference object.
        """
        try:
            pdb_id = reference[1]
            print(reference)
            self.logger.info(f"Retrieving PDB reference {pdb_id} from database")
            pdb_ref_exists = self.session.query(exists().where(PDBReference.pdb_id == pdb_id)).scalar()
            if not pdb_ref_exists:
                existing_sequence = self._get_or_create_sequence(sequence)
                return PDBReference(
                    pdb_id=pdb_id,
                    method=reference[2],
                    resolution=extract_float(reference[3]),
                    sequence=existing_sequence
                )
            else:
                return self.session.query(PDBReference).filter(PDBReference.pdb_id == pdb_id).first()
        except Exception as e:
            self.logger.error(f"Error retrieving or creating PDB reference for {pdb_id}: {e}")
            raise

    def _handle_go_reference(self, protein, reference):
        try:
            go_term = self._get_or_create_go_term(reference)
            association = self._get_or_create_association(protein.entry_name, go_term.go_id)
            if association is None:
                self.logger.info(
                    f"Association between {protein.entry_name} and GO Term {go_term.go_id} already exists.")
        except Exception as e:
            self.logger.error(f"Failed to handle GO reference for {protein.entry_name}: {e}")
            raise

    def _get_or_create_go_term(self, reference):
        """
        Retrieves or creates a Gene Ontology (GO) term in the database based on provided reference data.
        Args:
            reference (list): Contains the GO term details extracted from UniProt data.
        Returns:
            GOTerm: The retrieved or newly created GOTerm object.
        """
        try:
            go_id, category, description = reference[1], reference[2].split(':')[0], reference[2].split(':')[1]
            go_term = self.session.query(GOTerm).filter_by(go_id=go_id).first()
            if not go_term:
                go_term = GOTerm(go_id=go_id, category=category, description=description,
                                 evidences=reference[3].split(':')[0])
                self.session.add(go_term)
            return go_term
        except Exception as e:
            self.logger.error(f"Error retrieving or creating GO term {go_id}: {e}")
            raise

    def _get_or_create_association(self, entry_name, go_id):
        """
        Retrieves or creates an association between a protein and a Gene Ontology (GO) term in the database.
        It first checks if the association already exists to avoid duplication. If it does not exist, it creates
        a new association and adds it to the database.
        Args:
            entry_ptr (str): The name of the protein entry.
            go_id (str): The GO term identifier.
        Returns:
            ProteinGOTermAssociation: The newly created or existing association object.
        """
        try:
            exists_query = exists().where(
                ProteinGOTermAssociation.protein_entry_name == entry_name).where(
                ProteinGOTermAssociation.go_id == go_id)
            association_exists = self.session.query(exists_query).scalar()
            if not association_exists:
                association = ProteinGOTermAssociation(protein_entry_name=entry_name, go_id=go_id)
                self.session.add(association)
                return association
        except Exception as e:
            self.logger.error(f"Failed to retrieve or create association for {entry_name} and GO term {go_id}: {e}")
            raise