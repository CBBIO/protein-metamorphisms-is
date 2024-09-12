"""
PDB Extraction Tasks
=====================

The `PDBExtractor` class is responsible for managing the extraction, processing, and storage of protein structure data from the Protein Data Bank (PDB). It extends the `QueueTaskInitializer` class, leveraging RabbitMQ for task distribution and ensuring efficient handling of large-scale protein data.

**Purpose**

The `PDBExtractor` class manages the complete workflow of downloading, processing, and storing PDB data. This includes downloading PDB files, processing their models and chains, and saving the processed data into the system’s database for further analysis and research.

**Customization**

To create a custom PDB extraction task, subclass `PDBExtractor` and implement the `enqueue`, `process`, and `store_entry` methods. These methods define the logic for queuing tasks, processing them, and storing the results in the database.

**Key Features**

- **Directory Management**: Automatically creates directories for storing PDB files and processed models.
- **Task Enqueuing**: Dynamically enqueues tasks based on specific criteria like resolution thresholds.
- **Data Processing**: Downloads PDB files, processes each model and chain, and organizes results for storage.
- **Database Integration**: Ensures processed data is stored and organized in the database.
- **Error Handling**: Logs and handles errors during download and processing to ensure robustness.

**Example Usage**

Here is an example of how to subclass `PDBExtractor`:

.. code-block:: python

   from protein_metamorphisms_is.tasks.pdb import PDBExtractor

   class MyPDBExtractor(PDBExtractor):
       def enqueue(self):
           # Implementation of enqueue logic
           pass

       def process(self, pdb_id):
           # Implementation of process logic
           pass

       def store_entry(self, record):
           # Implementation of store_entry logic
           pass
"""

import gc
import os
import time
import traceback
import warnings

import gemmi
from Bio import PDB, SearchIO
from Bio.Data.PDBData import protein_letters_3to1_extended
from Bio.PDB import Select, MMCIFParser, MMCIFIO, FastMMCIFParser
from sqlalchemy.exc import NoResultFound
from sqlalchemy.sql.operators import or_

from protein_metamorphisms_is.tasks.queue import QueueTaskInitializer
from protein_metamorphisms_is.sql.model import PDBReference, Sequence, Structure, Model, PDBChains

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
        reference_attribute (str): The attribute used to reference PDB entries, typically 'PDB_ID'.
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
        self.reference_attribute = 'PDB_ID'
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
            if not os.path.exists(path):
                os.makedirs(path)
                self.logger.info(f"Created directory: {path}")
            else:
                self.logger.info(f"Directory already exists: {path}")

    def enqueue(self):
        """
        Enqueue tasks for processing based on PDB references.

        This method queries the database for PDB references that meet certain criteria,
        such as resolution threshold or method type, and enqueues them for processing.
        """
        resolution_threshold = self.conf.get("resolution_threshold")
        pdb_references = self.session.query(PDBReference).filter(
            or_(
                PDBReference.resolution <= resolution_threshold,
                PDBReference.method == "NMR"
            )
        ).all()
        for pdb_reference in pdb_references:
            self.logger.debug(f"Publishing task for accession code: {pdb_reference.pdb_id}")
            self.publish_task(pdb_reference.pdb_id)

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
            file_path = pdb_list.retrieve_pdb_file(pdb_id, file_format=self.conf.get("file_format", "mmCif"),
                                                   pdir=self.pdb_directory)
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

                        clean_file_path = os.path.join(self.models_directory,
                                                       f"{pdb_id}_clean_chain_{chain.name}_model_{model.name}.cif")
                        clean_structure.make_mmcif_document().write_file(clean_file_path)

                        data['chains'].append({
                            "chain_id": chain.name,
                            "sequence": ''.join(
                                [gemmi.find_tabulated_residue(residue.name).one_letter_code for residue in
                                 clean_chain]),
                            "model_id": model.name,
                            "file_path": clean_file_path
                        })

                        self.logger.info(f"Saved cleaned chain {chain.name} of model {model.name} to {clean_file_path}")

        except Exception as e:
            self.logger.error(f"Error processing PDB {pdb_id}: {e}\n{traceback.format_exc()}")
            return e

        return data

    def store_entry(self, record):
        """
        Store the processed PDB data into the database.

        This method saves the processed structures, models, and chains into the database,
        ensuring they are properly linked and retrievable.

        Args:
            record (dict): The processed data record to be stored.
        """
        try:
            pdb_structure = self.get_or_create_structure(record['pdb_id'])
            for chain_data in record['chains']:
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
            print(record)
            self.logger.error(f"Failed to store data in the database: {e}\n{traceback.format_exc()}")

    def get_or_create_structure(self, pdb_id, chain_id=None):
        """
        Retrieve or create a new structure in the database.

        This method checks if a structure already exists in the database; if not,
        it creates a new entry.

        Args:
            pdb_id (str): The PDB ID of the structure.
            chain_id (str, optional): The chain ID if applicable.

        Returns:
            Structure: The retrieved or newly created structure object.
        """
        file_name = f"{pdb_id}_chain_{chain_id}.cif" if chain_id else f"{pdb_id}.cif"

        existing_structure = self.session.query(Structure).filter_by(
            file_path=file_name).first()
        if not existing_structure:
            existing_structure = Structure(file_path=file_name)
            self.session.add(existing_structure)
            self.session.flush()

        return existing_structure

    def get_or_create_model(self, model_data, structure_id):
        """
        Retrieve or create a new model in the database.

        This method checks if a model already exists in the database; if not,
        it creates a new entry.

        Args:
            model_data (dict): Data of the model to be retrieved or created.
            structure_id (int): ID of the structure this model belongs to.

        Returns:
            Model: The retrieved or newly created model object.
        """
        model_id_str = model_data['model_id']
        file_name = os.path.basename(model_data['file_path'])

        existing_model = self.session.query(Model).filter_by(
            model_id=model_id_str, structure_id=structure_id
        ).first()
        if not existing_model:
            existing_model = Model(
                model_id=model_id_str,
                structure_id=structure_id,
                file_path=file_name
            )
            self.session.add(existing_model)
            self.session.flush()
        return existing_model

    def get_pdb_reference(self, pdb_id):
        """
        Retrieve the PDB reference from the database.

        Args:
            pdb_id (str): The PDB ID to be retrieved.

        Returns:
            PDBReference: The PDBReference object if found, else None.
        """
        try:
            pdb_reference = self.session.query(PDBReference).filter_by(pdb_id=pdb_id).one()
            return pdb_reference
        except NoResultFound:
            return None

    def get_or_create_chain(self, chain_data, pdb_reference_id, structure_id):
        """
        Retrieve or create a new chain in the database.

        This method checks if a chain already exists in the database; if not,
        it creates a new entry.

        Args:
            chain_data (dict): Data of the chain to be retrieved or created.
            pdb_reference_id (int): ID of the PDB reference.
            structure_id (int): ID of the structure this chain belongs to.

        Returns:
            PDBChains: The retrieved or newly created chain object.
        """
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
            existing_sequence = Sequence(sequence=sequence)
            self.session.add(existing_sequence)
            self.session.flush()
        return existing_sequence.id
