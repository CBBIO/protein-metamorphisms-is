import os

import mini3di
from Bio.PDB import MMCIFParser, Chain

from protein_metamorphisms_is.sql.model.entities.embedding.structure_3di import Structure3Di
from protein_metamorphisms_is.sql.model.entities.structure.state import State
from protein_metamorphisms_is.tasks.queue import QueueTaskInitializer


class Structure3DiManager(QueueTaskInitializer):
    """
    Manages the structural embedding process by encoding 3D atomic models into discrete 3Di representations.

    This class handles the loading and preprocessing of 3D structures from mmCIF files, applies the mini3di encoder
    to generate embeddings, and stores the resulting data in the database. The process is designed to run asynchronously
    within a distributed task queue system.

    Attributes
    ----------
    encoder : mini3di.Encoder
        Encoder object used to compute 3Di states from atomic structures.

    parser : Bio.PDB.MMCIFParser
        Parser for reading and interpreting mmCIF files.

    reference_attribute : str
        Name of the reference attribute used in task publication (default: 'model').
    """

    def __init__(self, conf):
        """
        Initializes the Structure3DiManager with the given configuration.

        This constructor sets up the encoder for 3Di state computation and the MMCIF parser
        for reading structural models. It inherits queue-based task management behavior from
        QueueTaskInitializer.

        Parameters
        ----------
        conf : dict
            Configuration dictionary. Must include:
            - 'data_directory': Path to the base directory where structural model files are stored.
            - 'limit_execution' (optional): Integer limit on the number of models to process (used for testing/debug).
        """

        super().__init__(conf)
        self.encoder = mini3di.Encoder()
        self.parser = MMCIFParser(QUIET=True)
        self.reference_attribute = "model"

    def enqueue(self):
        """
        Publishes structural embedding tasks for all available models (State entries).

        This method queries the database for all structural models stored in the `State` table.
        It then publishes a separate task for each model to be processed asynchronously.

        Behavior is optionally limited by the 'limit_execution' parameter in the configuration.

        Notes
        -----
        Each published task contains the dictionary representation of a `State` object,
        which includes model metadata such as file paths and identifiers.
        """

        states = self.session.query(State).all()
        if self.conf['limit_execution']:
            states = states[:self.conf['limit_execution']]
        for state in states:
            self.publish_task(state.__dict__)

    def process(self, model_info):
        """
        Loads and parses a 3D structural model, prepares it for encoding, and returns the embedding result.

        Given the metadata of a structural model (including its file path), this method loads
        the corresponding mmCIF file, extracts the first model available, and normalizes its
        residue chains. It then proceeds to encode the resulting chain using the 3Di encoder.

        Parameters
        ----------
        model_info : dict
            Dictionary containing information about the structural model, including:
            - 'file_path': Relative path to the mmCIF file.
            - 'id': Model identifier.

        Returns
        -------
        dict or None
            A dictionary with the model ID and its corresponding embedding, or None if an error occurs
            during parsing or encoding.

        Notes
        -----
        - Only the first model in the mmCIF file is processed.
        - If the structure file is empty or malformed, the method logs a warning and returns None.
        """

        file_path = os.path.join(self.conf['data_directory'], 'models', model_info['file_path'])
        try:
            structure = self.parser.get_structure("model", file_path)
            bio_model = next(structure.get_models())
        except StopIteration:
            self.logger.warning("No model found in the structure.")
            return

        new_chain = self.prepare_new_chain(bio_model)
        return self.process_chain(new_chain, model_info)

    def prepare_new_chain(self, bio_model):
        """
        Normalizes and rebuilds all chains in a structural model into a single renumbered chain.

        This method traverses all chains in the given model, copies and renumbers their residues
        sequentially, and merges them into a synthetic chain labeled 'A'. This ensures
        compatibility with the encoder, which expects a single continuous chain.

        Parameters
        ----------
        bio_model : Bio.PDB.Model.Model
            The structural model extracted from a parsed mmCIF file.

        Returns
        -------
        Bio.PDB.Chain.Chain
            A new chain object containing all residues from the original model, renumbered and unified.
        """

        residues_to_add = []
        for chain in list(bio_model.get_chains()):
            for residue in chain:
                new_id = (' ', len(residues_to_add) + 1, ' ')
                residue.id = new_id
                residues_to_add.append(residue.copy())
            bio_model.detach_child(chain.id)

        new_chain = Chain.Chain("A")
        for residue in residues_to_add:
            new_chain.add(residue)

        return new_chain

    def process_chain(self, chain, model_info):
        """
        Encodes a protein chain into a sequence of 3Di states and prepares the embedding result.

        This method applies the `mini3di.Encoder` to convert the given protein chain into a
        symbolic representation of its local structural patterns. If encoding fails, it logs
        detailed information about the chain to facilitate debugging.

        Parameters
        ----------
        chain : Bio.PDB.Chain.Chain
            Protein chain object prepared by `prepare_new_chain`.

        model_info : dict
            Dictionary with metadata about the original structure, including its 'id'.

        Returns
        -------
        dict or None
            Dictionary with keys:
            - 'model_id': identifier of the original structure.
            - 'embedding': encoded 3Di sequence.
            Returns None if encoding fails.
        """

        try:
            states = self.encoder.encode_chain(chain)
            sequence = self.encoder.build_sequence(states)
        except Exception as e:
            # Logging the entire chain might be too verbose, so log specific details instead
            self.logger.error(f"Error processing the chain: {e}", exc_info=True)
            return None

        embedding_result = {
            'model_id': model_info['id'],
            'embedding': sequence
        }
        return embedding_result

    def store_entry(self, record):
        """
        Persists the 3Di embedding of a structural model into the database.

        This method creates a new `Structure3Di` object using the model ID and the encoded
        embedding sequence. It attempts to store the entry in the database and handles
        any errors by rolling back the transaction and logging the exception.

        Parameters
        ----------
        record : dict
            Dictionary with the following keys:
            - 'model_id': Identifier of the structural model (`State.id`).
            - 'embedding': List or string representing the 3Di sequence.

        Notes
        -----
        - Commits the transaction only if insertion is successful.
        - Logs success or failure events accordingly.
        """

        new_embedding = Structure3Di(
            state_id=record['model_id'],
            embedding=record['embedding'],
        )
        try:
            self.session.add(new_embedding)
            self.session.commit()
            self.logger.info(f"Stored embedding for model ID {record['model_id']} successfully.")
        except Exception as e:
            self.session.rollback()
            self.logger.error(f"Failed to store embedding: {e}", exc_info=True)
