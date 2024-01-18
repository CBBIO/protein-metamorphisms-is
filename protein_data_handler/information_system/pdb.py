import traceback
from concurrent.futures import ThreadPoolExecutor

from Bio import PDB
from Bio.Data.PDBData import protein_letters_3to1
from Bio.PDB import Select, PDBIO, MMCIFParser
from Bio.PDB.Polypeptide import three_to_one
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql.operators import or_

from protein_data_handler.helpers.config.yaml import read_yaml_config
from protein_data_handler.information_system.base.bioinfo_extractor import BioinfoExtractorBase
from protein_data_handler.sql.model import PDBReference, PDBChains


class ChainSelect(Select):
    """
    A specialized selector for extracting specific chains from a PDB structure.
    This class is used in conjunction with Bio.PDB's PDBIO to selectively write
    specific chains of a PDB structure to a file.
    """

    def __init__(self, chain_id):
        """
        Initializes the ChainSelect with the specified chain ID

        Args:
            chain_id (str): The ID of the chain to be selected for output.
        """
        self.chain_id = chain_id

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
        return True


class PDBExtractor(BioinfoExtractorBase):
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
        Downloads and processes a given PDB structure using the Biopython library.

        This method handles the retrieval of PDB files from the PDB repository and
        processes them to extract relevant data, such as chain information and sequences.

        Args:
            pdb_reference (PDBReference): A PDBReference object containing PDB ID and other metadata.
        """

        local_session = sessionmaker(bind=self.engine)()

        pdb_id = pdb_reference.pdb_id
        pdbl = PDB.PDBList(server=self.conf.get("server", "ftp.wwpdb.org"), pdb=self.conf.get("pdb_path", "pdb_files"))
        try:
            file_path = pdbl.retrieve_pdb_file(pdb_id, file_format=self.conf.get("file_format", "pdb"),
                                               pdir=self.conf.get("pdb_path", "pdb_files"))
            self.logger.info(f"Descargado PDB {pdb_id}")

            # Procesar el archivo PDB descargado y poblar la base de datos
            self.populate_pdb_chains(file_path, pdb_reference.pdb_id, local_session)

        except Exception as e:
            self.logger.error(f"Error al descargar y procesar PDB {pdb_id}: {e}\n{traceback.format_exc()}")


    def populate_pdb_chains(self, pdb_file_path, pdb_reference_id, local_session):
        """
        Processes a PDB file to extract chain information and populates the database.

        This method parses the PDB file, extracts chain details and sequences, and
        saves them to the database. It also handles the creation of individual PDB
        files for each chain.

        Args:
            pdb_file_path (str): The file path to the PDB file.
            pdb_reference_id (int): The database ID for the PDB reference.
            local_session (Session): A SQLAlchemy session for database operations.
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
                sequence = ""
                for residue in chain:
                    if residue.id[0] == ' ' and residue.resname in protein_letters_3to1:  # Solo residuos estándar
                        sequence += protein_letters_3to1[residue.resname]
                pdb_chain = PDBChains(chains=chain_id, sequence=sequence, pdb_reference_id=pdb_reference_id_value)
                local_session.add(pdb_chain)

                # New code to write individual chain files
                chain_path = self.conf.get("pdb_chains_path", "pdb_chain_files")
                chain_file_path = f"{chain_path}/{pdb_reference_id}_{chain_id}.pdb"
                io = PDBIO()
                io.set_structure(structure)
                io.save(chain_file_path, select=ChainSelect(chain_id))

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
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            executor.map(self.download_and_process_pdb_structure, pdb_references)
