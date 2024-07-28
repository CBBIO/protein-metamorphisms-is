import os
import traceback
import warnings

from Bio import PDB, SearchIO
from Bio.Data.PDBData import protein_letters_3to1
from Bio.PDB import Select, MMCIFParser, MMCIFIO
from sqlalchemy.exc import NoResultFound
from sqlalchemy.sql.operators import or_

from protein_metamorphisms_is.base.queue import QueueTaskInitializer
from protein_metamorphisms_is.sql.model import PDBReference, Sequence, Structure, Model, PDBChains
import warnings
warnings.filterwarnings("ignore")




class ChainSelect(Select):
    def __init__(self, chain_id, model_id):
        self.chain_id = chain_id
        self.model_id = model_id

    def accept_chain(self, chain):
        return chain.get_id() == self.chain_id

    def accept_model(self, model):
        return model.get_id() == self.model_id


class PDBExtractor(QueueTaskInitializer):
    def __init__(self, conf):
        super().__init__(conf)
        self.reference_attribute = 'PDB_ID'
        self.data_directory = self.conf.get("data_directory", "/data")
        self.pdb_directory = os.path.join(self.data_directory, 'pdb')
        self.models_directory = os.path.join(self.data_directory, 'models')
        self.setup_directories()

    def setup_directories(self):
        for path in [self.pdb_directory, self.models_directory]:
            if not os.path.exists(path):
                os.makedirs(path)
                self.logger.info(f"Created directory: {path}")
            else:
                self.logger.info(f"Directory already exists: {path}")

    def enqueue(self):
        resolution_threshold = self.conf.get("resolution_threshold")
        pdb_references = self.session.query(PDBReference).filter(
            or_(
                PDBReference.resolution <= resolution_threshold,
                PDBReference.method == "NMR"
            )
        ).all()
        for pdb_reference in pdb_references:
            self.logger.debug(f"Publishing task for accession code: {pdb_reference.pdb_id}")
            print('pu')
            self.publish_task(pdb_reference.pdb_id)

    def process(self, pdb_id):
        pdb_list = PDB.PDBList(server=self.conf.get("server", "ftp.wwpdb.org"), pdb=self.pdb_directory)
        result = {
            "pdb_id": pdb_id,
            "models": [],
            "chains": []
        }

        try:
            file_path = pdb_list.retrieve_pdb_file(pdb_id, file_format=self.conf.get("file_format", "mmCif"),
                                                   pdir=self.pdb_directory)
            self.logger.info(f"Downloaded PDB {pdb_id} to {file_path}")
            parser = MMCIFParser()
            structure = parser.get_structure(pdb_id, file_path)
            io = MMCIFIO()

            for model in structure:
                model_file_path = os.path.join(self.models_directory, f"{pdb_id}_model_{model.id}.cif")
                io.set_structure(model)
                io.save(model_file_path)
                self.logger.info(f"Saved model {model.id} of PDB {pdb_id} to {model_file_path}")
                model_data = {
                    "model_id": model.id,
                    "file_path": model_file_path,
                    "chains": []
                }

                for chain in model:
                    chain_id = chain.get_id()
                    sequence = ""
                    for residue in chain:
                        if residue.id[0] == ' ' and residue.resname in protein_letters_3to1:
                            sequence += protein_letters_3to1[residue.resname]

                    chain_file_path = os.path.join(self.models_directory, f"{pdb_id}_chain_{chain_id}_model_{model.id}.cif")
                    if not os.path.isfile(chain_file_path):
                        io.set_structure(model)
                        io.save(chain_file_path, select=ChainSelect(chain_id, model.id))
                        self.logger.info(f"Saved chain {chain_id} of model {model.id} to {chain_file_path}")
                    model_data['chains'].append({
                        "chain_id": chain_id,
                        "sequence": sequence,
                        "model_id": model.id,
                        "file_path": chain_file_path
                    })

                result['models'].append(model_data)

        except Exception as e:
            self.logger.error(f"Error downloading and processing PDB {pdb_id}: {e}\n{traceback.format_exc()}")
            return e

        return result

    def store_entry(self, record):
        print(record)
        try:
            pdb_structure = self.get_or_create_structure(record['pdb_id'])

            for model_data in record['models']:
                model = self.get_or_create_model(model_data, pdb_structure.id)

                for chain_data in model_data['chains']:
                    chain_structure = self.get_or_create_structure(record['pdb_id'], chain_data['chain_id'])
                    model = self.get_or_create_model(chain_data, chain_structure.id)
                    pdb_reference = self.get_pdb_reference(record['pdb_id'])

                    if pdb_reference:
                        pdb_reference.structure_id = pdb_structure.id
                        self.session.add(pdb_reference)

                    self.get_or_create_chain(chain_data, pdb_reference.id, chain_structure.id)

            self.session.commit()
        except Exception as e:
            self.session.rollback()
            self.logger.error(f"Failed to store data in the database: {e}\n{traceback.format_exc()}")
        finally:
            self.session.close()

    def get_or_create_structure(self, pdb_id, chain_id=None):
        file_name = f"{pdb_id}_chain_{chain_id}.cif" if chain_id else f"{pdb_id}.cif"

        existing_structure = self.session.query(Structure).filter_by(
            file_path=file_name).first()  # Guarda solo el nombre del archivo
        if not existing_structure:
            existing_structure = Structure(file_path=file_name)  # Guarda solo el nombre del archivo
            self.session.add(existing_structure)
            self.session.flush()

        return existing_structure

    def get_or_create_model(self, model_data, structure_id):
        model_id_str = str(model_data['model_id'])
        file_name = os.path.basename(model_data['file_path'])

        existing_model = self.session.query(Model).filter_by(
            model_id=model_id_str, structure_id=structure_id
        ).first()
        if not existing_model:
            existing_model = Model(
                model_id=model_id_str,
                structure_id=structure_id,
                file_path=file_name  # Guarda solo el nombre del archivo
            )
            self.session.add(existing_model)
            self.session.flush()
        return existing_model

    def get_pdb_reference(self, pdb_id):
        try:
            pdb_reference = self.session.query(PDBReference).filter_by(pdb_id=pdb_id).one()
            return pdb_reference
        except NoResultFound:
            return None

    def get_or_create_chain(self, chain_data, pdb_reference_id, structure_id):
        sequence_id = self.get_or_create_sequence(chain_data['sequence'])
        existing_chain = self.session.query(PDBChains).filter_by(
            chains=chain_data['chain_id'], pdb_reference_id=pdb_reference_id, structure_id=structure_id
        ).first()
        if not existing_chain:
            new_chain = PDBChains(
                chains=chain_data['chain_id'],
                sequence_id=sequence_id,
                pdb_reference_id=pdb_reference_id,
                structure_id=structure_id
            )
            self.session.add(new_chain)
            self.session.flush()
        return existing_chain

    def get_or_create_sequence(self, sequence):
        existing_sequence = self.session.query(Sequence).filter_by(sequence=sequence).first()
        if not existing_sequence:
            existing_sequence = Sequence(sequence=sequence)
            self.session.add(existing_sequence)
            self.session.flush()
        return existing_sequence.id
