import traceback
from http.client import HTTPException

from Bio import ExPASy, SwissProt
from sqlalchemy import exists

from protein_metamorphisms_is.sql.model.entities.go_annotation.go_annotation import ProteinGOTermAnnotation
from protein_metamorphisms_is.sql.model.entities.go_annotation.go_term import GOTerm
from protein_metamorphisms_is.sql.model.entities.protein.accesion import Accession
from protein_metamorphisms_is.sql.model.entities.protein.protein import Protein
from protein_metamorphisms_is.sql.model.entities.sequence.sequence import Sequence
from protein_metamorphisms_is.sql.model.entities.structure.structure import Structure
from protein_metamorphisms_is.tasks.queue import QueueTaskInitializer
from protein_metamorphisms_is.helpers.parser.parser import extract_float


class UniProtExtractor(QueueTaskInitializer):
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
        self.reference_attribute = 'Accession'

    def enqueue(self):
        accessions = self.session.query(Accession).all()
        for accession in accessions:
            self.logger.debug(f"Publishing task for {self.reference_attribute} code: {accession.code}")
            self.publish_task(accession.code)

    # Existing method in uniprot.py

    def process(self, accession_code):
        """
        Download detailed protein information from UniProt using ExPASy and SwissProt.
        ExPASy is a Bioinformatics Resource Portal which provides access to scientific databases and software tools,
        while SwissProt is a manually annotated and reviewed protein sequence database part of UniProt.
        Args:
            accession_dict (dict): The dictionary containing the accession code of the protein record to download.
        Returns:
            record: A dictionary containing detailed protein information.
        Raises:
            ValueError: If no SwissProt record is found for the accession code.
        """
        if not accession_code:
            raise ValueError("No accession code provided in the input dictionary.")

        try:
            handle = ExPASy.get_sprot_raw(accession_code)
            record = SwissProt.read(handle)
            return record
        except ValueError:
            error_message = f"No SwissProt record found for accession code {accession_code}"
            self.logger.error(error_message)
            raise ValueError(error_message)
        except Exception as e:
            error_message = f"Error downloading the entry {accession_code}: {e}"
            self.logger.error(error_message)
            raise Exception(error_message)

    def store_entry(self, data):
        """
        Store the UniProt data in the database, ensuring that Accession is linked to the correct Protein.
        """
        try:
            self.logger.debug(f"Storing data for entry_name {data.entry_name}.")

            # Retrieve or create the protein
            protein = self.get_or_create_protein(data)
            self.handle_cross_references(protein, data.cross_references)

            # Update the details of the protein
            self.update_protein_details(protein, data)

            # Link Accession to Protein
            accession_entry = (
                self.session.query(Accession)
                .filter_by(code=data.accessions[0])  # Assuming `data.accessions[0]` is the primary accession
                .first()
            )

            if not accession_entry:
                accession_entry = Accession(
                    code=data.accessions[0],
                    protein_id=protein.id,  # Link the protein ID
                    primary=True
                )
                self.session.add(accession_entry)
                self.logger.debug(f"Created new accession {data.accessions[0]} linked to protein {protein.id}.")
            else:
                # Update the existing accession if necessary
                if accession_entry.protein_id != protein.id:
                    accession_entry.protein_id = protein.id
                    self.logger.debug(f"Updated accession {data.accessions[0]} with new protein ID {protein.id}.")

            # Commit changes
            self.session.commit()
            self.logger.info(f"Successfully stored accession {data.accessions[0]} for protein {protein.id}.")

        except Exception as e:
            self.session.rollback()
            self.logger.error(f"Failed to store data for entry {data.entry_name}: {e}")
            raise

    def get_or_create_protein(self, data):
        """
        Retrieves an existing protein record from the database using the entry_name from data or creates a new one if it does not exist.
        Args:
            data (SwissProt.Record): Data containing the entry_name of the protein.
        Returns:
            Protein: The retrieved or newly created protein object.
        """
        try:
            protein = self.session.query(Protein).filter_by(id=data.entry_name).first()
            if not protein:
                protein = Protein(id=data.entry_name)
                self.session.add(protein)
                self.logger.debug(f"Created new protein record for entry_name {data.entry_name}.")
            else:
                self.logger.debug(f"Retrieved existing protein record for entry_name {data.entry_name}.")
            return protein
        except Exception as e:
            self.logger.error(
                f"Failed to retrieve or create protein for entry_name {data.entry_name}: {e}\n{traceback.format_exc()}")
            raise e

    def update_protein_details(self, protein, data):
        """
        Updates the protein details in the database and ensures the protein is linked to the correct sequence.

        Args:
            protein (Protein): The protein object to update.
            data (SwissProt.Record): The data containing the new protein information.
        """
        try:
            # Get or create the associated sequence
            sequence = self.get_or_create_sequence(data.sequence)

            # Flush the session to ensure the sequence.id is available
            self.session.flush()

            # Update the protein's sequence link if different
            if protein.sequence_id != sequence.id:
                protein.sequence_id = sequence.id

            # Update other protein fields
            protein.data_class = data.data_class
            protein.molecule_type = data.molecule_type
            protein.sequence_length = data.sequence_length
            protein.created_date = data.created[0] if data.created else None
            protein.sequence_update_date = data.sequence_update[0] if data.sequence_update else None
            protein.annotation_update_date = data.annotation_update[0] if data.annotation_update else None
            protein.description = data.description
            protein.gene_name = str(data.gene_name)
            protein.organism = data.organism
            protein.organelle = data.organelle
            protein.organism_classification = ",".join(
                data.organism_classification) if data.organism_classification else None
            protein.taxonomy_id = ",".join(data.taxonomy_id) if data.taxonomy_id else None
            protein.host_organism = ",".join(data.host_organism) if data.host_organism else None
            protein.host_taxonomy_id = ",".join(data.host_taxonomy_id) if data.host_taxonomy_id else None
            protein.comments = "; ".join(data.comments) if data.comments else None
            protein.keywords = data.keywords
            protein.protein_existence = data.protein_existence
            protein.seqinfo = data.seqinfo

            # Ensure the protein object is added to the session for update
            self.session.add(protein)
        except Exception as e:
            self.logger.error(f"Failed to update protein details: {e}")
            raise e

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
                self.logger.debug(f"Created new sequence record.")
            else:
                self.logger.debug(f"Retrieved existing sequence record.")
            return existing_sequence
        except Exception as e:
            self.logger.error(f"Failed to retrieve or create sequence: {e}\n{traceback.format_exc()}")
            raise e

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
            self.logger.error(
                f"Failed to handle cross references for {protein.id}: {e}\n{traceback.format_exc()}")
            raise e

    def handle_pdb_reference(self, protein, reference):
        """
        Specific handler for PDB references, extracting relevant segment if specified and storing it.
        Args:
            protein (Protein): The protein object associated with the PDB reference.
            reference (list): PDB reference data, including PDB ID, method, and resolution.
        """
        try:
            if reference[0] == "PDB":
                structure = self.get_or_create_structure(reference, protein.id)
                self.session.add(structure)
                self.logger.debug(f"Linked structure {reference[1]} to protein {protein.id}.")
        except Exception as e:
            self.logger.error(f"Error handling PDB reference {reference}: {e}\n{traceback.format_exc()}")
            raise e

    def get_or_create_structure(self, reference, protein_id):
        """
        Retrieves or creates a Structure entry in the database based on the provided PDB ID.
        If the structure does not exist, it creates a new one using the associated reference data.
        Args:
            reference (list): PDB reference data including PDB ID, method, and resolution.
            protein_id (str): The ID of the associated protein.
        Returns:
            Structure: The retrieved or newly created Structure object.
        """
        try:
            pdb_id = reference[1]
            method = reference[2]
            resolution = extract_float(reference[3])

            self.logger.info(f"Retrieving structure {pdb_id} from the database.")
            structure_exists = self.session.query(exists().where(Structure.id == pdb_id)).scalar()

            if not structure_exists:
                self.logger.debug(f"Creating new structure record for {pdb_id}.")
                structure = Structure(
                    id=pdb_id,
                    protein_id=protein_id,
                    method=method,
                    resolution=resolution,
                    file_path=f"structures/{pdb_id}.cif"  # Adjust the path as needed
                )
                return structure
            else:
                self.logger.debug(f"Retrieved existing structure record for {pdb_id}.")
                return self.session.query(Structure).filter(Structure.id == pdb_id).first()
        except Exception as e:
            self.logger.error(f"Error retrieving or creating structure for {pdb_id}: {e}\n{traceback.format_exc()}")
            raise e

    def handle_go_reference(self, protein, reference):
        try:
            evidence = reference[3].split(":")[0] if reference[
                3] else "UNKNOWN"  # Usa un valor por defecto si falta evidence
            if evidence not in self.conf['allowed_evidences']:
                self.logger.debug(
                    f"Skipping GO reference {reference[1]} for {protein.id} due to disallowed evidence type.")
                return

            # Crear o recuperar el término GO
            go_term = self.get_or_create_go_term(reference)

            # Verificar si la asociación ya existe
            association = self.get_or_create_association(protein.id, go_term.go_id, evidence_code=evidence)

            if association is None:
                self.logger.info(f"Association between {protein.id} and GO Term {go_term.go_id} already exists.")
            else:
                self.logger.debug(f"Created new association between {protein.id} and GO Term {go_term.go_id}.")
        except Exception as e:
            self.logger.error(f"Failed to handle GO reference for {protein.id}: {e}\n{traceback.format_exc()}")
            raise e

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
                self.logger.debug(f"Creating new GO term record for {go_id}.")
                go_term = GOTerm(go_id=go_id, category=category, description=description)
                self.session.add(go_term)
            else:
                self.logger.debug(f"Retrieved existing GO term record for {go_id}.")
            return go_term
        except Exception as e:
            self.logger.error(f"Error retrieving or creating GO term {go_id}: {e}\n{traceback.format_exc()}")
            raise e

    def get_or_create_association(self, entry_name, go_id, evidence_code, is_transferred=False,
                                  source_cluster_id=None, target_cluster_id=None, distance=None,
                                  embedding_type_id=None):
        """
        Create or retrieve a GOAnnotation entry.
        """
        try:
            # Verificar si la asociación ya existe
            exists_query = exists().where(
                ProteinGOTermAnnotation.protein_id == entry_name).where(
                ProteinGOTermAnnotation.go_id == go_id)

            association_exists = self.session.query(exists_query).scalar()

            if not association_exists:
                self.logger.debug(f"Creating new association for {entry_name} and GO term {go_id}.")
                # Crear la asociación con el evidence_code proporcionado
                association = ProteinGOTermAnnotation(
                    protein_id=entry_name,
                    go_id=go_id,
                    evidence_code=evidence_code  # Agregar evidence_code
                )
                self.session.add(association)
                return association
            else:
                self.logger.debug(f"Association between {entry_name} and GO term {go_id} already exists.")
                return None
        except Exception as e:
            self.logger.error(
                f"Failed to retrieve or create association for {entry_name} and GO term {go_id}: {e}")
            raise e
