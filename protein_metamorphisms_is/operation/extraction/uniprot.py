import traceback

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
    The UniProtExtractor class is responsible for extracting, processing, and storing protein data
    from the UniProt database. It extends the QueueTaskInitializer to handle task queuing and
    database integration seamlessly.

    **Purpose**

    The UniProtExtractor class provides functionalities to manage and process protein data
    from UniProt, ensuring it is properly stored and cross-referenced in the database.

    **Key Features**

    - **Task Management**: Supports enqueuing tasks for processing accession codes from UniProt.
    - **Protein Data Retrieval**: Downloads detailed protein information using BioPython's ExPASy and SwissProt modules.
    - **Database Integration**: Ensures extracted data is stored and linked correctly in the database.
    - **Cross-Reference Handling**: Processes references to PDB structures and Gene Ontology (GO) terms.

    **Example Usage**

    Below is an example demonstrating how to use the `UniProtExtractor` class:

    .. code-block:: python

        from protein_metamorphisms_is.tasks.uniprot import UniProtExtractor

        # Configuration for UniProtExtractor
        config = {
            'allowed_evidences': ['EXP', 'IDA'],
                }

        # Initialize and start the extractor
        uniprot_extractor = UniProtExtractor(config)
        uniprot_extractor.start()

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
        """
        Reads all accession codes from the database and publishes them as tasks to be processed.
        """
        try:
            # Recuperar accesiones de la base de datos
            accessions = self.session.query(Accession).all()

            if not accessions:
                self.logger.info("No accession codes found in the database. Enqueue operation skipped.")
                return

            # Aplicar el límite de ejecución si está configurado
            limit_execution = self.conf.get("limit_execution")
            if limit_execution and isinstance(limit_execution, int):
                accessions = accessions[:limit_execution]

            self.logger.info(f"Retrieved {len(accessions)} accession codes from the database.")

            # Publicar tareas
            for accession in accessions:
                self.logger.debug(f"Publishing task for accession code '{accession.code}'...")
                self.publish_task(accession.code)

            self.logger.info(f"Successfully enqueued {len(accessions)} tasks for accession codes.")
        except Exception as e:
            self.logger.error(
                f"Error occurred during enqueue operation: {e}\n{traceback.format_exc()}"
            )
            raise

    def process(self, accession_code):
        """
        Downloads detailed protein information from UniProt using ExPASy and SwissProt modules.

        Args:
            accession_code (str): The accession code of the protein record to download.

        Returns:
            SwissProt.Record: The protein record fetched from UniProt.

        Raises:
            ValueError: If no SwissProt record is found for the accession code.
            Exception: If any other error occurs during data retrieval.
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
        except Exception:
            error_message = f"Error downloading UniProt entry for accession '{accession_code}'"
            self.logger.error(error_message)
            raise Exception(error_message)

    def store_entry(self, data):
        """
        Stores the retrieved UniProt data into the database.

        Links accession codes to the associated protein and ensures all references and annotations are updated.

        Args:
            data (SwissProt.Record): The UniProt data to store.

        Raises:
            Exception: If the storage process fails.
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
            self.logger.info(
                f"Accession '{data.accessions[0]}' successfully linked to protein '{protein.id}'."
            )

        except Exception as e:
            self.session.rollback()
            self.logger.error(f"Failed to store data for entry {data.entry_name}: {e}")
            raise

    def get_or_create_protein(self, data):
        """
        Retrieves or creates a protein record in the database.

        Args:
            data (SwissProt.Record): The UniProt record containing protein details.

        Returns:
            Protein: The retrieved or newly created protein object.

        Raises:
            Exception: If an error occurs during the operation.
        """
        try:
            protein = self.session.query(Protein).filter_by(id=data.entry_name).first()
            if not protein:
                protein = Protein(id=data.entry_name)
                self.session.add(protein)
                self.logger.debug(
                    f"New protein record created for UniProt entry '{data.entry_name}'."
                )

            else:
                self.logger.debug(f"Retrieved existing protein record for entry_name {data.entry_name}.")
            return protein
        except Exception as e:
            self.logger.error(
                f"Failed to retrieve or create protein for entry_name {data.entry_name}")
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
                self.logger.debug("Created new sequence record.")
            else:
                self.logger.debug("Retrieved existing sequence record.")
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
            go_id, category, description, _ = reference[1], reference[2].split(":")[0], reference[2].split(":")[
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
