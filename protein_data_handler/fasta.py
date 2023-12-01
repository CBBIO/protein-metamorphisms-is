import os
import requests
from concurrent.futures import ThreadPoolExecutor
import logging
from requests.exceptions import RequestException

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


class FastaHandler:
    """
    Clase para descargar archivos FASTA de la base de datos de PDB
        (Protein Data Bank).

    Esta clase permite descargar archivos FASTA, que contienen secuencias de
        aminoácidos o nucleótidos, para un conjunto de identificadores de
        PDB proporcionados.

    :param session: Sesión de requests utilizada para realizar las descargas.
    :type session: requests.Session
    :param data_dir: Directorio donde se guardarán los archivos FASTA
        descargados.
    :type data_dir: str
    """

    def __init__(self, session, data_dir, output_dir):
        """
        Inicializa el descargador de FASTA con una sesión de base de datos y un
            directorio de datos.

        :param session: Sesión de requests para realizar las descargas.
        :param data_dir: Ruta del directorio donde se almacenarán los archivos
            FASTA.
        """
        self.session = session
        self.data_dir = data_dir
        self.output_dir = output_dir

        if not os.path.exists(data_dir):
            os.makedirs(data_dir, exist_ok=True)

        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)

    def download_fastas(self, pdb_ids, max_workers=10):
        """
        Descarga archivos FASTA para un conjunto de IDs de PDB utilizando
            múltiples hilos.

        :param pdb_ids: Lista de identificadores de PDB para los cuales
            descargar los archivos FASTA.
        :type pdb_ids: list[str]
        :param max_workers: Número máximo de hilos para usar en la descarga.
        :type max_workers: int
        :raises ValueError: Si `pdb_ids` no es una lista de cadenas de texto.
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

        :param pdb_id: Identificador de PDB para el cual descargar el archivo
            FASTA.
        :type pdb_id: str
        :raises ValueError: Si `pdb_id` no es una cadena de texto.
        :raises RequestException: Si ocurre un error en la solicitud HTTP.
        :raises IOError: Si ocurre un error al escribir el archivo descargado.
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

    def merge_fastas(self, pdb_ids, merge_name):
        """
        Combina archivos FASTA para un conjunto de IDs de PDB en un único
        archivo. Descarga los archivos FASTA que no están presentes en el
        directorio local.

        :param pdb_ids: Lista de identificadores de PDB para los cuales
            combinar los archivos FASTA.
        :type pdb_ids: list[str]
        """
        missing_files = []
        for pdb_id in pdb_ids:
            file_path = os.path.join(self.data_dir, f"{pdb_id}.fasta")
            if not os.path.isfile(file_path):
                missing_files.append(pdb_id)

        if missing_files:
            logging.info("Descargando archivos FASTA faltantes.")
            self.download_fastas(missing_files)

        with (open(os.path.join(self.output_dir, f"{merge_name}.fasta"), 'w')
              as outfile):
            for pdb_id in pdb_ids:
                file_path = os.path.join(self.data_dir, f"{pdb_id}.fasta")
                if os.path.isfile(file_path):
                    with open(file_path, 'r') as infile:
                        outfile.write(infile.read())
                        outfile.write('\n')
                else:
                    logging.warning(f"Archivo no encontrado: {file_path}")
