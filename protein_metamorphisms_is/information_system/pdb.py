import os
import traceback
from concurrent.futures import ThreadPoolExecutor

from Bio import PDB
from Bio.Data.PDBData import protein_letters_3to1
from Bio.PDB import Select, MMCIFParser, MMCIFIO
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql.operators import or_

from protein_metamorphisms_is.information_system.base.extractor import ExtractorBase
from protein_metamorphisms_is.sql.model import PDBReference, PDBChains


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

    def start(self):
        """
        Begins the process of extracting data from the PDB.

        This method initiates the download and processing of PDB structures based on
        predefined criteria (like resolution threshold) from the configuration.
        """
        try:
            self.logger.info("Iniciando la extracción de datos del PDB")

            pdb_references = self.load_pdb_ids()

            self.download_pdb_structures(pdb_references)

        except Exception as e:
            self.logger.error(f"Error durante el proceso de extracción: {e}")

    def load_pdb_ids(self):
        """
        Load PDB IDs from the database that meet the specified resolution threshold.

        This method queries the database for PDB entries with resolution values
        below a certain threshold, indicating higher quality structures.

        Returns:
            list: A list of PDBReference objects meeting the resolution criteria.
        """
        resolution_threshold = self.conf.get("resolution_threshold")  # Corregido aquí
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

        Steps:
            1. Initiates a database session.
            2. Retrieves the PDB file based on the PDB ID from the PDBReference object.
            3. Downloads the file to a specified directory in the desired format.
            4. Calls 'populate_pdb_chains' to process the downloaded file and store chain information
               in the database.
        """
        local_session = sessionmaker(bind=self.engine)()

        pdb_id = pdb_reference.pdb_id
        pdbl = PDB.PDBList(server=self.conf.get("server", "ftp.wwpdb.org"), pdb=self.conf.get("pdb_path", "pdb_files"))
        try:
            file_path = pdbl.retrieve_pdb_file(pdb_id, file_format=self.conf.get("file_format", "mmCif"),
                                               pdir=self.conf.get("pdb_path", "pdb_files"))
            self.logger.info(f"Descargado PDB {pdb_id}")

            self.populate_pdb_chains(file_path, pdb_reference.pdb_id, local_session)

        except Exception as e:
            self.logger.error(f"Error al descargar y procesar PDB {pdb_id}: {e}\n{traceback.format_exc()}")

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

        The method proceeds as follows:
        - It queries the database to check if the provided 'pdb_reference_id' exists.
        - For each chain in the PDB file:
            - Extracts the chain ID, model ID, and the amino acid sequence.
            - Creates a new PDBChains object with the extracted data and adds it to the session.
            - Generates a CIF file for the chain, storing it in a specified directory if it doesn't already exist.

        Note:
        - The method assumes 'protein_letters_3to1' is a dictionary mapping three-letter amino acid
          codes to their single-letter counterparts.
        - 'ChainSelect' is a custom selector used by MMCIFIO for saving individual chains.
        - The directory for saving individual chain CIF files is configurable and defaults to 'pdb_chain_files'.
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
                pdb_chain = PDBChains(chains=chain_id, sequence=sequence, pdb_reference_id=pdb_reference_id_value,
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

    def download_pdb_structures(self, pdb_references):
        """
        Downloads and processes PDB structures in parallel given their IDs.

        This method uses concurrent processing to handle multiple PDB structure
        downloads and processing simultaneously, enhancing efficiency.

        Args:
            pdb_references (list): A list of PDBReference objects to be processed.
        """
        max_workers = self.conf.get("max_workers", 5)
        self.logger.info("Downloading PDB structures with {} workers".format(max_workers))
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            executor.map(self.download_and_process_pdb_structure, pdb_references)
