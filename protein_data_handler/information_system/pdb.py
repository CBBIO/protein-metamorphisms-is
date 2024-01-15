from concurrent.futures import ThreadPoolExecutor

from Bio import PDB
from Bio.PDB import Select, PDBIO, MMCIFParser
from sqlalchemy.orm import sessionmaker

from protein_data_handler.helpers.config.yaml import read_yaml_config
from protein_data_handler.information_system.base.bioinfo_extractor import BioinfoExtractorBase
from protein_data_handler.sql.model import PDBReference, PDBChains


class ChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.get_id() == self.chain_id

    def accept_model(self, model):
        # Aceptar solo el modelo con el ID especificado
        return True


class PDBExtractor(BioinfoExtractorBase):
    """
    A class for extracting data from UniProt.

    This class extends BioinfoExtractorBase and provides specific implementations
    for extracting and processing data from UniProt.

    Args:
        conf (dict): Configuration dictionary containing necessary parameters.
    """

    def __init__(self, conf):
        """
        Initialize the UniProtExtractor class.

        Sets up the configuration and logger. Initializes the database session if required.
        """
        super().__init__(conf, session_required=True)

    def start(self):
        """
        Comienza el proceso de extracción de datos del PDB.
        """
        try:
            self.logger.info("Iniciando la extracción de datos del PDB")

            pdb_references = self.load_pdb_ids()

            self.download_pdb_structures(pdb_references)


        except Exception as e:
            self.logger.error(f"Error durante el proceso de extracción: {e}")

    def load_pdb_ids(self):
        """
        Load PDB IDs from the database that meet the resolution threshold.
        """
        resolution_threshold = self.conf.get("resolution_threshold")  # Corregido aquí
        pdb_references = self.session.query(PDBReference).filter(
            PDBReference.resolution <= resolution_threshold).all()
        return pdb_references

    def download_and_process_pdb_structure(self, pdb_reference):
        """
        Descarga y procesa una estructura PDB dada su ID utilizando Biopython.
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
            self.logger.error(f"Error al descargar y procesar PDB {pdb_id}: {e}")

    def populate_pdb_chains(self, pdb_file_path, pdb_reference_id,local_session):
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
                sequence = "".join([residue.get_resname() for residue in chain if residue.id[0] == ' '])
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
        Descarga y procesa estructuras PDB en paralelo dados sus IDs.
        """
        max_workers = self.conf.get("max_workers", 5)
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            executor.map(self.download_and_process_pdb_structure, pdb_references)
