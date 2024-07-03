import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from http.client import HTTPException

import pandas as pd
import requests
from Bio import ExPASy, SwissProt
from rq import Queue
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

    def fetch(self):
        """
        Initiates the download of UniProt entries using Redis queues to manage the download tasks.
        """
        self.logger.info("Starting the download of UniProt entries.")

        accessions = self.session.query(Accession).all()
        self.logger.info(f"Total proteins to download: {len(accessions)}")

        max_workers = self.conf.get("max_workers", 10)  # Default to 10 if not specified in config

        for accession in accessions:
            job = self.process_queue.enqueue(download_record, accession.accession_code)






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
            protein = self.get_or_create_protein(data)
            self.update_protein_details(protein, data)
            self.handle_accessions(protein, data.accessions)
            self.handle_cross_references(protein, data.cross_references)
            self.session.commit()
            self.logger.info(f"Successfully stored or updated data for protein {data.entry_name}.")
        except HTTPException as e:
            self.logger.error(f"HTTP error occurred while processing {data.entry_name}: {e}")
            self.session.rollback()
        except Exception as e:
            self.logger.error(f"An unexpected error occurred while processing {data.entry_name}: {e}")
            self.session.rollback()

    def get_or_create_protein(self, data):
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

    def update_protein_details(self, protein, data):
        """
        Updates the protein details based on the provided SwissProt data. All related protein attributes like sequence,
        molecule type, organism details, etc., are updated.
        Args:
            protein (Protein): The protein object to update.
            data (SwissProt.Record): Data used to update the protein object.
        """
        protein.sequence = self.get_or_create_sequence(data.sequence)
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

    def get_or_create_sequence(self, sequence):
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

    def handle_accessions(self, protein, accessions):
        """
        Processes and links accession codes to the specified protein, ensuring no duplicates in the database.
        Args:
            protein (Protein): The protein object to which the accession codes will be linked.
            accessions (list): List of accession codes.
        """
        try:
            for accession_code in accessions:
                accession = self.get_or_create_accession(accession_code)
                accession.protein_entry_name = protein.entry_name
                self.session.add(accession)
        except Exception as e:
            self.logger.error(f"Error handling accession {accession_code} from {accessions}: {e}")
            raise

    def get_or_create_accession(self, accession_code):
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

    def handle_cross_references(self, protein, cross_references):
        """
        Manages the cross-references associated with the protein, such as database links to PDB and GO terms.
        Args:
            protein (Protein): The protein object to manage.
            cross_references (list): List of cross-reference data.
        """
        try:
            for reference in cross_references:
                if reference[0] == "PDB":
                    self.handle_pdb_reference(protein, reference)
                elif reference[0] == "GO":
                    self.handle_go_reference(protein, reference)
        except Exception as e:
            self.logger.error(f"Failed to handle cross references for {protein.entry_name}: {e}")
            raise

    def handle_pdb_reference(self, protein, reference):
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
                        pdb_ref = self.get_or_create_pdb_reference(reference, sequence)
                        pdb_ref.protein = protein
                        self.session.add(pdb_ref)
        except Exception as e:
            self.logger.error(f"Error handling PDB reference {reference}: {e}")
            raise

    def get_or_create_pdb_reference(self, reference, sequence):
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
            self.logger.info(f"Retrieving PDB reference {pdb_id} from database")
            pdb_ref_exists = self.session.query(exists().where(PDBReference.pdb_id == pdb_id)).scalar()
            if not pdb_ref_exists:
                existing_sequence = self.get_or_create_sequence(sequence)
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

    def handle_go_reference(self, protein, reference):
        try:
            evidence = reference[3].split(":")[0]
            if evidence not in self.conf['allowed_evidences']:
                return
            go_term = self.get_or_create_go_term(reference)
            association = self.get_or_create_association(protein.entry_name, go_term.go_id)
            if association is None:
                self.logger.info(
                    f"Association between {protein.entry_name} and GO Term {go_term.go_id} already exists.")
        except Exception as e:
            self.logger.error(f"Failed to handle GO reference for {protein.entry_name}: {e}")
            raise

    def get_or_create_go_term(self, reference):
        """
        Retrieves or creates a Gene Ontology (GO) term in the database based on provided reference data.
        Args:
            reference (list): Contains the GO term details extracted from UniProt data.
        Returns:
            GOTerm: The retrieved or newly created GOTerm object.
        """
        try:
            go_id, category, description, evidence = reference[1], reference[2].split(":")[0], reference[2].split(":")[
                1], reference[3].split(":")[0]
            go_term = self.session.query(GOTerm).filter_by(go_id=go_id).first()

            if not go_term:
                go_term = GOTerm(go_id=go_id, category=category, description=description,
                                 evidences=reference[3].split(':')[0])
                self.session.add(go_term)
            return go_term
        except Exception as e:
            self.logger.error(f"Error retrieving or creating GO term {go_id}: {e}")
            raise

    def get_or_create_association(self, entry_name, go_id):
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


def download_record(accession_code):
    """
    Download detailed protein information from UniProt using ExPASy and SwissProt.
    ExPASy is a Bioinformatics Resource Portal which provides access to scientific databases and software tools,
    while SwissProt is a manually annotated and reviewed protein sequence database part of UniProt.
    Args:
        accession_code (str): The accession code of the protein record to download.
    Returns:
        record: A SwissProt record object containing detailed protein information.
    Raises:
        Exception: Raises an exception with a message indicating what went wrong during the download.
    """
    try:
        handle = ExPASy.get_sprot_raw(accession_code)
        record = SwissProt.read(handle)
        return record
    except Exception as e:
        raise Exception(f"Error downloading the entry {accession_code}: {e}")