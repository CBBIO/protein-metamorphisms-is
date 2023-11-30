import os
import requests
from concurrent.futures import ThreadPoolExecutor
import logging
from requests.exceptions import RequestException

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


class FastaDownloader:
    """
    Clase para descargar archivos FASTA de la base de datos de PDB.
    """

    def __init__(self, session, data_dir):
        """
        Inicializa el descargador de FASTA con una sesión de base de datos.
        """
        self.session = session
        self.data_dir = data_dir
        if not os.path.exists(data_dir):
            os.makedirs(data_dir, exist_ok=True)

    def download_fastas(self, pdb_ids, max_workers=10):
        """
        Descarga archivos FASTA para un conjunto de IDs de PDB
        utilizando múltiples hilos.
        """
        logging.info(f"Descarga de {len(pdb_ids)} estructuras FASTA.")
        if not isinstance(pdb_ids, list) or not all(isinstance(id, str)
                                                    for id in pdb_ids):
            raise ValueError("pdb_ids debe ser una lista de cadenas de texto.")

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            executor.map(self.download_fasta, pdb_ids)

    def download_fasta(self, pdb_id):
        """
        Descarga un archivo FASTA individual de la base de datos de PDB.
        """

        if not isinstance(pdb_id, str):
            raise ValueError("pdb_id debe ser una cadena de texto.")

        url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
        try:
            response = requests.get(url)
            response.raise_for_status()

            file_path = os.path.join(self.data_dir, f"{pdb_id}.fasta")
            with open(file_path, "w") as file:
                file.write(response.text)
            logging.info(f"FASTA descargado para {pdb_id} en {file_path}")

        except RequestException as e:
            logging.error(f"Error al descargar FASTA para {pdb_id}: {e}")
        except IOError as e:
            logging.error(f"Error al escribir el archivo para {pdb_id}: {e}")
