import os
import traceback
from concurrent.futures import ThreadPoolExecutor

from Bio import PDB
from Bio.Data.PDBData import protein_letters_3to1
from Bio.PDB import Select, MMCIFParser, MMCIFIO
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql.operators import or_

from protein_metamorphisms_is.information_system.base.extractor import ExtractorBase
from protein_metamorphisms_is.sql.model import PDBReference, PDBChains, Sequence


class ChainSelect(Select):
    """
    A specialized selector for extracting specific chains from a PDB structure.
    This class is used in conjunction with Bio.PDB's PDBIO to selectively write
    specific chains of a PDB structure to a file.
    """

    def __init__(self, chain_id, model_id):
        """
        Initializes the ChainSelect with the specified chain ID

        Args:
            chain_id (str): The ID of the chain to be selected for output.
        """
        self.chain_id = chain_id
        self.model_id = model_id

    def accept_chain(self, chain):
        """
        Determines if a chain should be accepted (written to file).

        Args:
            chain (Chain): A chain object from a PDB structure.

        Returns:
            bool: True if the chain's ID matches the desired chain_id, False otherwise.
        """
        return chain.get_id() == self.chain_id

    def accept_model(self, model):
        """
        Accepts all models. In PDB files, models represent different conformations
        of the structure, commonly used in NMR structures

        Args:
            model (Model): A model object from a PDB structure

        Returns:
            bool: Always True, as all models are accepted.
        """
        return model.get_id() == self.model_id


class PDBExtractor(ExtractorBase):
    """
    A class for extracting and processing PDB (Protein Data Bank) structures.

    This class extends BioinfoExtractorBase, providing specific implementations
    for downloading, parsing, and processing structural data from the PDB,
    a global repository of information about the 3D structures of large biological
    molecules, including proteins and nucleic acids.

    Args:
        conf (dict): Configuration dictionary containing necessary parameters.
    """

    def __init__(self, conf):
        """
        Initialize the PDBExtractor class.

        Sets up the configuration and logger. Initializes the database session if required.
        """
        super().__init__(conf, session_required=True)

    def set_targets(self):
        """
        Sets the targets for extraction, which are the PDB IDs.
        """
        self.logger.info("Loading PDB IDs from the database")
        self.pdb_references = self.load_pdb_ids()
        self.logger.info(f"Loaded {len(self.pdb_references)} PDB IDs")

    def fetch(self):
        """
        Fetches the data by downloading PDB structures in parallel.
        """
        max_workers = self.conf.get("max_workers", 5)
        self.logger.info(f"Downloading PDB structures with {max_workers} workers")
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            executor.map(self.download_and_process_pdb_structure, self.pdb_references)

    def store_entry(self, data):
        """
        Stores a single entry of processed data into the database.

        Args:
            data (tuple): A tuple containing (pdb_file_path, pdb_reference_id).
        """
        pdb_file_path, pdb_reference_id = data
        local_session = sessionmaker(bind=self.engine)()
        self.populate_pdb_chains(pdb_file_path, pdb_reference_id, local_session)

    def load_pdb_ids(self):
        """
        Load PDB IDs from the database that meet the specified resolution threshold.

        This method queries the database for PDB entries with resolution values
        below a certain threshold, indicating higher quality structures.

        Returns:
            list: A list of PDBReference objects meeting the resolution criteria.
        """
        resolution_threshold = self.conf.get("resolution_threshold")
        pdb_references = self.session.query(PDBReference).filter(
            or_(
                PDBReference.resolution <= resolution_threshold,
                PDBReference.method == "NMR"
            )
        ).all()
        return pdb_references

    def download_and_process_pdb_structure(self, pdb_reference):
        """
        Downloads and processes a PDB structure using the Biopython library.

        This method is responsible for downloading PDB files from the specified PDB repository,
        then processing these files to extract relevant data. It primarily focuses on retrieving
        chain information and sequences from the PDB structure. The process involves two main steps:
        downloading the PDB file and then populating the database with chain details extracted from
        the file.

        Args:
            pdb_reference (PDBReference): A PDBReference object that contains the PDB ID and other
                                          related metadata necessary for downloading and processing the file.
        """
        pdb_id = pdb_reference.pdb_id
        pdbl = PDB.PDBList(server=self.conf.get("server", "ftp.wwpdb.org"), pdb=self.conf.get("pdb_path", "pdb_files"))
        try:
            file_path = pdbl.retrieve_pdb_file(pdb_id, file_format=self.conf.get("file_format", "mmCif"),
                                               pdir=self.conf.get("pdb_path", "pdb_files"))
            self.logger.info(f"Downloaded PDB {pdb_id}")
            self.data_queue.put((file_path, pdb_reference.pdb_id))
        except Exception as e:
            self.logger.error(f"Error downloading and processing PDB {pdb_id}: {e}\n{traceback.format_exc()}")

    def populate_pdb_chains(self, pdb_file_path, pdb_reference_id, local_session):
        """
        Processes a PDB file, extracting and storing chain information in the database.

        This method uses the MMCIFParser to parse the PDB file specified by 'pdb_file_path'.
        It then extracts details such as chain identifiers, models, and sequences from the file.
        Each chain's information, along with its associated sequence and reference to the PDB file,
        is stored in the database. Additionally, this method generates individual CIF files for
        each chain in the structure, which are saved to a specified directory.

        Args:
            pdb_file_path (str): Path to the PDB file to be processed.
            pdb_reference_id (int): The unique database identifier for the PDB reference.
            local_session (Session): An active SQLAlchemy session for executing database operations.
        """
        parser = MMCIFParser()
        structure = parser.get_structure(pdb_reference_id, pdb_file_path)
        pdb_reference_id_result = local_session.query(PDBReference.id).filter(
            PDBReference.pdb_id == pdb_reference_id).first()

        if pdb_reference_id_result is not None:
            (pdb_reference_id_value,) = pdb_reference_id_result
        else:
            pdb_reference_id_value = None

        for model in structure:
            for chain in model:
                chain_id = chain.get_id()
                model_id = model.get_id()
                sequence = ""
                for residue in chain:
                    if residue.id[0] == ' ' and residue.resname in protein_letters_3to1:
                        sequence += protein_letters_3to1[residue.resname]
                existing_sequence = local_session.query(Sequence).filter_by(sequence=sequence).first()
                if not existing_sequence:
                    existing_sequence = Sequence(sequence=sequence)
                    local_session.add(existing_sequence)
                pdb_chain = PDBChains(chains=chain_id, sequence=existing_sequence,
                                      pdb_reference_id=pdb_reference_id_value,
                                      model=model_id)
                local_session.add(pdb_chain)

                chain_path = self.conf.get("pdb_chains_path", "pdb_chain_files")
                chain_file_path = f"{chain_path}/{pdb_reference_id}_{chain_id}_{model.get_id()}.cif"

                if not os.path.isfile(chain_file_path):
                    io = MMCIFIO()
                    io.set_structure(structure)
                    io.save(chain_file_path, select=ChainSelect(chain_id, model_id))
        local_session.commit()
        local_session.close()
