import os
import traceback
import warnings

import gemmi
from Bio import PDB
from Bio.PDB import Select
from sqlalchemy import func, text
from sqlalchemy.exc import NoResultFound
from sqlalchemy.sql import or_

from protein_metamorphisms_is.sql.model.entities.protein.protein import Protein
from protein_metamorphisms_is.sql.model.entities.sequence.sequence import Sequence
from protein_metamorphisms_is.tasks.queue import QueueTaskInitializer
from protein_metamorphisms_is.sql.model.entities.structure.structure import Structure
from protein_metamorphisms_is.sql.model.entities.structure.chain import Chain
from protein_metamorphisms_is.sql.model.entities.structure.state import State

warnings.filterwarnings("ignore")


class ChainSelect(Select):
    """
    A helper class used during PDB file parsing to filter chains and models.
    """

    def __init__(self, chain_id, model_id):
        self.chain_id = chain_id
        self.model_id = model_id

    def accept_chain(self, chain):
        """
        Filter to accept only the specified chain.
        """
        return chain.get_id() == self.chain_id

    def accept_model(self, model):
        """
        Filter to accept only the specified model.
        """
        return model.get_id() == self.model_id


class PDBExtractor(QueueTaskInitializer):
    """
    The PDBExtractor class manages the extraction and processing of protein structures
    from the Protein Data Bank (PDB). It handles the entire workflow from downloading
    PDB files, processing each model and chain, to storing the processed data in the
    system’s database.

    Attributes:
        reference_attribute (str): The attribute used to reference PDB entries, typically 'id'.
        data_directory (str): The root directory where data files are stored. Defaults to '/data'.
        pdb_directory (str): Directory where raw PDB files are stored.
        models_directory (str): Directory where processed model files are saved.
    """

    def __init__(self, conf):
        """
        Initialize the PDBExtractor.

        This constructor sets up the necessary directories for storing PDB files and
        processed models and initializes the configuration.

        Args:
            conf (dict): Configuration dictionary containing settings like data directories.
        """
        super().__init__(conf)
        self.reference_attribute = 'id'  # Actualmente 'id' en la clase Structure
        self.data_directory = self.conf.get("data_directory", "/data")
        self.pdb_directory = os.path.join(self.data_directory, 'pdb')
        self.models_directory = os.path.join(self.data_directory, 'models')
        self.setup_directories()

    def setup_directories(self):
        """
        Set up the directories needed for PDB file storage and processing.

        Creates the directories for storing raw PDB files and processed model files
        if they do not already exist.
        """
        for path in [self.pdb_directory, self.models_directory]:
            try:
                if not os.path.exists(path):
                    os.makedirs(path)
                    self.logger.info(f"Created directory: {path}")
                else:
                    self.logger.info(f"Directory already exists: {path}")
            except OSError as e:
                self.logger.error(f"Failed to create directory {path}: {e}")
                raise

    def enqueue(self):
        """
        Enqueue tasks for processing based on Structure entries.

        This method queries the database for Structure entries that meet certain criteria,
        such as resolution threshold or method type, and enqueues them for processing.
        """
        resolution_threshold = self.conf.get("resolution_threshold")
        if resolution_threshold is None:
            self.logger.error("Resolution threshold not defined in configuration.")
            return

        try:
            structures = self.session.query(Structure).filter(
                or_(
                    Structure.resolution <= resolution_threshold,
                    Structure.method == "NMR"
                )
            ).all()
            for structure in structures:
                self.logger.debug(f"Publishing task for PDB ID: {structure.id}")
                self.publish_task(structure.id)
        except Exception as e:
            self.logger.error(f"Error enqueuing tasks: {e}\n{traceback.format_exc()}")

    def process(self, pdb_id):
        """
        Process a PDB file, extracting and saving its models and chains.

        This method downloads the PDB file, processes each model and chain,
        and saves the results in a structured format.

        Args:
            pdb_id (str): The PDB ID to be processed.

        Returns:
            dict: A dictionary containing the processed data.
        """
        pdb_list = PDB.PDBList(server=self.conf.get("server", "ftp.wwpdb.org"), pdb=self.pdb_directory)
        data = {
            "pdb_id": pdb_id,
            "chains": []
        }

        # Lista de aminoácidos a conservar
        amino_acids = [
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS',
            'GLU', 'GLY', 'HIS', 'ILE', 'LEU',
            'LYS', 'MET', 'PHE', 'PRO', 'SER',
            'THR', 'TRP', 'TYR', 'VAL'
        ]

        try:
            file_path = pdb_list.retrieve_pdb_file(
                pdb_id,
                file_format=self.conf.get("file_format", "mmCif"),
                pdir=self.pdb_directory
            )
            self.logger.info(f"Downloaded PDB {pdb_id} to {file_path}")

            structure = gemmi.read_structure(file_path)
            structure.remove_ligands_and_waters()
            structure.remove_hydrogens()
            structure.remove_empty_chains()

            for model in structure:
                for chain in model:
                    # Crear una nueva estructura limpia
                    clean_structure = gemmi.Structure()
                    clean_model = gemmi.Model(model.name)
                    clean_chain = gemmi.Chain(chain.name)

                    # Iterar sobre los residuos y filtrar los que sean aminoácidos
                    for residue in chain:
                        if residue.name in amino_acids:
                            clean_chain.add_residue(residue)

                    # Si la cadena tiene residuos después de filtrar, guárdala
                    if len(clean_chain) > 0:
                        clean_model.add_chain(clean_chain)
                        clean_structure.add_model(clean_model)

                        clean_file_path = os.path.join(
                            self.models_directory,
                            f"{pdb_id}_clean_chain_{chain.name}_model_{model.name}.cif"
                        )
                        clean_structure.make_mmcif_document().write_file(clean_file_path)

                        # Generar secuencia de una letra
                        sequence = ''.join(
                            [gemmi.find_tabulated_residue(residue.name).one_letter_code for residue in clean_chain]
                        )

                        data['chains'].append({
                            "chain_id": chain.name,
                            "sequence": sequence,
                            "file_path": clean_file_path,
                            'model': model.name
                        })

                        self.logger.info(f"Saved cleaned chain {chain.name} of model {model.name} to {clean_file_path}")

        except Exception as e:
            self.logger.error(f"Error processing PDB {pdb_id}: {e}\n{traceback.format_exc()}")
            return {"error": str(e), "traceback": traceback.format_exc()}

        return data

    def store_entry(self, record):
        """
        Store the processed PDB data into the database.

        This method saves the processed states and chains into the database,
        ensuring they are properly linked to an existing structure.

        Args:
            record (dict): The processed data record to be stored.
        """
        pdb_id = record["pdb_id"]
        chains_data = record["chains"]

        try:
            # Propagate the existing Structure
            structure = self.session.query(Structure).filter_by(id=pdb_id).one()

            # Iterate over chains to store Chain, Sequence, and State
            for chain_data in chains_data:
                chain_id = chain_data["chain_id"]
                sequence = chain_data["sequence"]
                file_path = chain_data["file_path"]
                # Handle the Sequence
                sequence_entry = self.session.query(Sequence).filter_by(sequence=sequence).one_or_none()
                if sequence_entry is None:
                    # Create new Sequence entry if it doesn't exist
                    sequence_entry = Sequence(sequence=sequence)
                    self.session.add(sequence_entry)

                # Handle the Chain
                chain = self.session.query(Chain).filter_by(name=chain_id, structure_id=structure.id).one_or_none()
                if chain is None:
                    # Create new Chain entry if it doesn't exist
                    file_path = chain_data["file_path"]

                    chain = Chain(name=chain_id, structure_id=structure.id, sequence_id=sequence_entry.id)
                    self.session.add(chain)




                # Handle the State
                state = self.session.query(State).filter_by(chain_id=chain.id, structure_id=structure.id, model_id=chain_data['model']).one_or_none()
                if state is None:
                    state = State(
                        model_id=chain_data['model'],  # Use model from the chain data
                        file_path=file_path,
                        chain_id=chain.id,
                        structure_id=structure.id
                    )
                    self.session.add(state)

            # Commit all changes to the database
            self.session.commit()
            self.logger.info(f"Stored data for PDB ID: {pdb_id}")

        except Exception as e:
            self.session.rollback()
            self.logger.error(f"Failed to store data for PDB ID: {pdb_id}: {e}")
            raise
