import os
import traceback
import warnings

import gemmi
from Bio import PDB
from Bio.PDB import Select
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
                            "file_path": clean_file_path
                        })

                        self.logger.info(f"Saved cleaned chain {chain.name} of model {model.name} to {clean_file_path}")

        except Exception as e:
            self.logger.error(f"Error processing PDB {pdb_id}: {e}\n{traceback.format_exc()}")
            return {"error": str(e), "traceback": traceback.format_exc()}

        return data

    def store_entry(self, record):
        """
        Store the processed PDB data into the database.

        This method saves the processed structures, chains, and states into the database,
        ensuring they are properly linked and retrievable.

        Args:
            record (dict): The processed data record to be stored.
        """
        if 'error' in record:
            self.logger.error(f"Record processing failed: {record['error']}")
            return

        try:
            # Obtener o crear la estructura
            structure = self.get_or_create_structure(record['pdb_id'])
            for chain_data in record['chains']:
                # Obtener o crear la cadena
                chain = self.get_or_create_chain(record['pdb_id'], chain_data['chain_id'])

                # Obtener o crear la secuencia
                sequence = self.get_or_create_sequence(chain_data['sequence'])

                # Obtener o crear el estado (State) en lugar del modelo
                state = self.get_or_create_state(chain_data['file_path'], chain.id, structure.id)

                # Si deseas asociar la secuencia con el estado o cadena, asegúrate de hacerlo aquí
                # Por ejemplo, si la cadena tiene una relación con la secuencia:
                # chain.sequence_id = sequence.id
                # Y si el State tiene una relación con la secuencia:
                # state.sequence_id = sequence.id

            self.session.commit()

        except Exception as e:
            self.session.rollback()
            self.logger.error(f"Failed to store data in the database: {e}\n{traceback.format_exc()}")
            self.logger.debug(f"Record data: {record}")

    def get_or_create_structure(self, pdb_id):
        """
        Retrieve or create a new structure in the database.

        This method checks if a structure already exists in the database; if not,
        it creates a new entry.

        Args:
            pdb_id (str): The PDB ID of the structure.

        Returns:
            Structure: The retrieved or newly created structure object.
        """
        try:
            structure = self.session.query(Structure).filter_by(id=pdb_id).one()
            return structure
        except NoResultFound:
            # Asumiendo que 'protein_id' debe ser asignado de alguna manera
            # Necesitas obtener o definir 'protein_id' según tu lógica
            protein_id = self.get_protein_id(pdb_id)  # Implementa este método según corresponda
            if not protein_id:
                raise ValueError(f"No Protein found for PDB ID {pdb_id}")

            new_structure = Structure(
                id=pdb_id,
                protein_id=protein_id,
                method=self.conf.get("default_method", "X-ray"),  # Asigna valores predeterminados o extrae de los datos
                resolution=self.conf.get("default_resolution", 2.0),  # Asigna valores predeterminados o extrae de los datos
                file_path=f"{pdb_id}.cif"
            )
            self.session.add(new_structure)
            self.session.flush()
            return new_structure

    def get_or_create_chain(self, pdb_id, chain_id):
        """
        Retrieve or create a new chain in the database.

        This method checks if a chain already exists in the database; if not,
        it creates a new entry.

        Args:
            pdb_id (str): The PDB ID of the structure.
            chain_id (str): The chain ID to be retrieved or created.

        Returns:
            Chain: The retrieved or newly created chain object.
        """
        try:
            chain = self.session.query(Chain).filter_by(
                id=chain_id,
                structure_id=pdb_id
            ).one()
            return chain
        except NoResultFound:
            new_chain = Chain(
                id=chain_id,
                structure_id=pdb_id
            )
            self.session.add(new_chain)
            self.session.flush()
            return new_chain

    def get_or_create_sequence(self, sequence):
        """
        Retrieve or create a new sequence in the database.

        This method checks if a sequence already exists in the database; if not,
        it creates a new entry.

        Args:
            sequence (str): The sequence to be retrieved or created.

        Returns:
            int: The ID of the retrieved or newly created sequence.
        """
        existing_sequence = self.session.query(Sequence).filter_by(sequence=sequence).first()
        if not existing_sequence:
            new_sequence = Sequence(sequence=sequence)
            self.session.add(new_sequence)
            self.session.flush()
            return new_sequence.id
        return existing_sequence.id

    def get_or_create_state(self, file_path, chain_id, structure_id):
        """
        Retrieve or create a new state in the database.

        This method checks if a state already exists in the database based on the file_path,
        chain_id, and structure_id; if not, it creates a new entry.

        Args:
            file_path (str): The file path of the processed chain.
            chain_id (str): The chain ID associated with the state.
            structure_id (str): The structure ID associated with the state.

        Returns:
            State: The retrieved or newly created state object.
        """
        try:
            state = self.session.query(State).filter_by(
                file_path=file_path,
                chain_id=chain_id,
                structure_id=structure_id
            ).one()
            return state
        except NoResultFound:
            new_state = State(
                file_path=file_path,
                chain_id=chain_id,
                structure_id=structure_id
            )
            self.session.add(new_state)
            self.session.flush()
            return new_state

    def get_protein_id(self, pdb_id):
        """
        Retrieve the protein_id associated with a given PDB ID.

        Implementa este método según la lógica de tu aplicación para obtener
        el 'protein_id' correspondiente al 'pdb_id'.

        Args:
            pdb_id (str): The PDB ID.

        Returns:
            str: The associated protein_id.
        """
        # Ejemplo de implementación:
        # Esto depende de cómo estés gestionando las relaciones entre proteínas y estructuras.
        protein = self.session.query(Protein).filter_by(pdb_id=pdb_id).first()
        if protein:
            return protein.id
        return None
