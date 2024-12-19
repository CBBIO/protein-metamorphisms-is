import os

import mini3di
from Bio.PDB import MMCIFParser, Chain

from protein_metamorphisms_is.sql.model.entities.embedding.structure_3di import Structure3Di
from protein_metamorphisms_is.sql.model.entities.structure.state import State
from protein_metamorphisms_is.tasks.queue import QueueTaskInitializer


class Structure3DiManager(QueueTaskInitializer):
    def __init__(self, conf):
        super().__init__(conf)
        self.encoder = mini3di.Encoder()
        self.parser = MMCIFParser(QUIET=True)
        self.reference_attribute = "model"

    def enqueue(self):
        states = self.session.query(State).all()
        if self.conf['limit_execution']:
            states = states[:self.conf['limit_execution']]
        for state in states:
            self.publish_task(state.__dict__)

    def process(self, model_info):
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
        try:
            states = self.encoder.encode_chain(chain)
            sequence = self.encoder.build_sequence(states)
        except Exception as e:
            # Logging the entire chain might be too verbose, so log specific details instead
            self.log_chain_details(chain, model_info)
            self.logger.error(f"Error processing the chain: {e}", exc_info=True)
            return None

        embedding_result = {
            'model_id': model_info['id'],
            'embedding': sequence
        }
        return embedding_result

    def log_chain_details(self, chain, model_info):
        try:
            print(model_info)
            print(f"Chain ID: {chain.id}")
            print(f"Number of residues: {len(list(chain.get_residues()))}")

            # Print each residue's details along with all atom coordinates and names
            for residue in chain.get_residues():
                print(f"Residue ID: {residue.id}, Residue Name: {residue.resname}")
                for atom in residue.get_atoms():
                    print(f"Atom Name: {atom.name}, Coordinates: {atom.coord}")

        except Exception as log_error:
            self.logger.error(f"Failed during logging chain details: {log_error}", exc_info=True)

    def store_entry(self, record):
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
