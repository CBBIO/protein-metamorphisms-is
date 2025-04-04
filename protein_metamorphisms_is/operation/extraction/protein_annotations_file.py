import requests
from protein_metamorphisms_is.sql.model.entities.go_annotation.go_annotation import ProteinGOTermAnnotation
from protein_metamorphisms_is.sql.model.entities.go_annotation.go_term import GOTerm
from protein_metamorphisms_is.sql.model.entities.protein.protein import Protein
from protein_metamorphisms_is.sql.model.entities.sequence.sequence import Sequence
from protein_metamorphisms_is.tasks.queue import QueueTaskInitializer


class GOAnnotationsQueueProcessor(QueueTaskInitializer):
    def __init__(self, conf):
        super().__init__(conf)
        self.file_path = self.conf['goa_annotations_file']

    def enqueue(self):
        """
        Enqueue tasks from the GO annotations file, limited by `limit_execution` if configured.
        This will send protein entry IDs and GO terms to the processing queue.
        """
        self.logger.info(f"Enqueueing tasks from file: {self.file_path}")
        try:
            with open(self.file_path, 'r') as f:
                lines = f.readlines()

            # Apply the limit if it's specified in the configuration
            limit_execution = self.conf.get("limit_execution")
            if limit_execution and isinstance(limit_execution, int):
                lines = lines[:limit_execution]
                self.logger.info(f"Limiting to the first {limit_execution} entries.")

            for line in lines:
                line = line.strip()  # Remove any leading/trailing whitespace
                if line:
                    # Extract protein entry ID and GO terms
                    parts = line.split('\t')
                    if len(parts) == 2:
                        protein_entry_id = parts[0]
                        go_terms = parts[1].split(',')

                        # Enqueue each task for processing
                        task_data = {
                            'protein_entry_id': protein_entry_id,
                            'go_terms': go_terms
                        }
                        self.publish_task(task_data)
                        self.logger.info(f"Task enqueued for protein {protein_entry_id}.")

        except Exception as e:
            self.logger.error(f"Error enqueuing tasks: {e}")
            raise

    def process(self, data):
        """
        Process a single task: get sequence and process GO terms for a given protein.
        """
        protein_entry_id = data['protein_entry_id']
        go_terms = data['go_terms']

        sequence = self.get_sequence_from_external_source(protein_entry_id)

        # Return a structured object containing protein data and annotations
        result = {
            'protein': protein_entry_id,
            'go_terms': go_terms,
            'sequence': sequence
        }

        # Return the structured data to be processed by store_entry
        return result

    def store_entry(self, data):
        """
        Stores the protein data, sequence, and GO annotations into the database.

        Args:
            data (dict): Contains protein entry ID, sequence, and associated GO terms.

        Raises:
            Exception: If the storage process fails.
        """
        try:
            # Log the start of the storing process
            self.logger.debug(f"Storing data for protein entry {data['protein']}.")

            # Step 1: Retrieve or create the protein entry
            protein = self.get_or_create_protein(data['protein'])

            # Step 2: Create or retrieve the sequence for the protein
            sequence = self.get_or_create_sequence(data['sequence'])

            # Link sequence to protein if not already linked
            if protein.sequence_id != sequence.id:
                protein.sequence_id = sequence.id
                self.logger.debug(f"Linked sequence to protein {protein.id}.")

            # Step 3: Process GO terms and store them
            for go_term in data['go_terms']:
                go_term_entry = self.get_or_create_go_term(go_term)

                # Create or retrieve the Protein-GO Term annotation
                self.get_or_create_association(protein.id, go_term_entry.go_id)

            # Commit all changes to the database
            self.session.commit()
            self.logger.info(f"Protein {protein.id} successfully updated with sequence and GO terms.")

        except Exception as e:
            # If anything fails, rollback the changes
            self.session.rollback()
            self.logger.error(f"Failed to store data for protein entry {data['protein']}: {e}")
            raise

    def get_or_create_sequence(self, sequence):
        """
        Retrieves or creates a sequence entity in the database.

        Args:
            sequence (str): The amino acid sequence of the protein.

        Returns:
            Sequence: The retrieved or newly created sequence object.
        """
        try:
            # Check if the sequence already exists in the database
            existing_sequence = self.session.query(Sequence).filter_by(sequence=sequence).first()
            if not existing_sequence:
                # Create a new sequence record if it doesn't exist
                existing_sequence = Sequence(sequence=sequence)
                self.session.add(existing_sequence)
                self.logger.debug(
                    f"Created new sequence record for sequence {sequence[:10]}...")  # Log a portion of the sequence for brevity
            else:
                self.logger.debug(f"Found existing sequence record for sequence {sequence[:10]}...")

            return existing_sequence
        except Exception as e:
            self.logger.error(f"Error retrieving or creating sequence: {e}")
            raise e

    def get_or_create_association(self, protein_id, go_id, evidence_code="UNKNOWN"):
        """
        Create or retrieve a Protein-GO Term annotation in the database.

        Args:
            protein_id (int): The ID of the protein.
            go_id (str): The GO term ID.
            evidence_code (str): The evidence code for the GO annotation. Defaults to "UNKNOWN".

        Returns:
            ProteinGOTermAnnotation: The retrieved or newly created association.
        """
        try:
            # Check if the association already exists
            existing_association = self.session.query(ProteinGOTermAnnotation).filter_by(protein_id=protein_id,
                                                                                         go_id=go_id).first()
            if not existing_association:
                # Create a new association record if it doesn't exist
                association = ProteinGOTermAnnotation(protein_id=protein_id, go_id=go_id, evidence_code=evidence_code)
                self.session.add(association)
                self.logger.debug(
                    f"Created new association for protein {protein_id} and GO term {go_id} with evidence code {evidence_code}.")
            else:
                self.logger.debug(f"Association already exists for protein {protein_id} and GO term {go_id}.")

            return existing_association
        except Exception as e:
            self.logger.error(f"Error retrieving or creating GO term association: {e}")
            raise e

    def get_sequence_from_external_source(self, protein_entry_id):
        """
        Retrieve the protein sequence from UniProt using the protein entry ID.
        """
        url = f"https://rest.uniprot.org/uniprotkb/{protein_entry_id}.json"
        try:
            response = requests.get(url)
            response.raise_for_status()  # Raise error for bad HTTP responses
            data = response.json()  # Assuming you have the response data

            # Extract the sequence safely
            sequence = data.get('sequence', {}).get('value', '')  # Extract the sequence

            if not sequence:
                self.logger.error(f"No sequence found for protein entry ID {protein_entry_id}")
                return None

            return sequence
        except requests.exceptions.RequestException as e:
            self.logger.error(f"Failed to fetch sequence for {protein_entry_id} from UniProt: {e}")
            return None

    def get_or_create_protein(self, protein_entry_id):
        """
        Retrieve or create a protein record in the database using the protein entry ID.
        """
        protein = self.session.query(Protein).filter_by(id=protein_entry_id).first()
        if not protein:
            # Create a new protein record if it doesn't exist
            protein = Protein(id=protein_entry_id)
            self.session.add(protein)
            self.session.commit()
            self.logger.debug(f"Created new protein record for protein entry ID {protein_entry_id}.")
        else:
            self.logger.debug(f"Protein record found for protein entry ID {protein_entry_id}.")

        return protein

    def get_or_create_go_term(self, go_term):
        """
        Retrieve or create a GO term entry in the database.
        """
        go_term_entry = self.session.query(GOTerm).filter_by(go_id=go_term).first()
        if not go_term_entry:
            go_term_entry = GOTerm(go_id=go_term)
            self.session.add(go_term_entry)
            self.session.commit()
            self.logger.debug(f"Created new GO term entry for GO term {go_term}.")
        return go_term_entry
