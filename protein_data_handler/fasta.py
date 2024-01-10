import os
import requests
from concurrent.futures import ThreadPoolExecutor
import logging

from requests.exceptions import RequestException

from protein_data_handler.helpers.parser.parser import extract_and_parse_fasta
from protein_data_handler.models.uniprot import PDBChains, PDBReference, Cluster

from pycdhit import cd_hit, read_clstr

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
        """
        logging.info(f"Descarga de {len(pdb_ids)} estructuras FASTA.")

        if not isinstance(pdb_ids, list) or not all(isinstance(id, str) for id in pdb_ids):
            raise ValueError("pdb_ids debe ser una lista de cadenas de texto.")

        to_download = []
        already_downloaded = []

        for pdb_id in pdb_ids:
            file_path = os.path.join(self.data_dir, f"{pdb_id}.fasta")
            if os.path.exists(file_path):
                already_downloaded.append(pdb_id)
            else:
                to_download.append(pdb_id)

        logging.info(
            f"Ya en disco: {len(already_downloaded)} archivos. Necesitan descarga: {len(to_download)} archivos.")

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            results = executor.map(self.download_and_store_fasta, pdb_ids)
        for file_path, chains in results:
            if file_path and chains:
                for chain in chains:
                    pdb_id, chain_number, chain_id, sequence = chain
                    pdb_reference_id = self.session.query(PDBReference).filter_by(pdb_id=pdb_id).first().id

                    # Verificar si la cadena ya existe en la base de datos
                    existing_chain = self.session.query(PDBChains).filter_by(
                        pdb_reference_id=pdb_reference_id,
                        chain_number=chain_number,
                        chains=chain_id
                    ).first()

                    if not existing_chain:
                        pdb_chain = PDBChains(pdb_reference_id=pdb_reference_id, chain_number=chain_number,
                                              chains=chain_id,
                                              sequence=sequence)
                        self.session.add(pdb_chain)

        self.session.commit()

    def download_and_store_fasta(self, pdb_id):
        """
        Descarga un archivo FASTA individual de la base de datos de PDB.
        Retorna el path del archivo y los datos extraídos si la descarga fue exitosa.
        """
        if not isinstance(pdb_id, str):
            raise ValueError("pdb_id debe ser una cadena de texto.")

        file_path = os.path.join(self.data_dir, f"{pdb_id}.fasta")
        if not os.path.exists(file_path):
            url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
            try:
                response = requests.get(url)
                response.raise_for_status()

                with open(file_path, "w") as file:
                    file.write(response.text)
                logging.info(f"FASTA descargado para {pdb_id} en {file_path}")

            except RequestException as e:
                logging.error(f"Error al descargar FASTA para {pdb_id}: {e}")
                return None, None
            except IOError as e:
                logging.error(f"Error al escribir el archivo para {pdb_id}: {e}")
                return None, None

        chains = extract_and_parse_fasta(file_path)
        return file_path, chains

    def cluster_fastas(self, input_file, threshold=0.7):
        """
        Agrupa archivos FASTA utilizando cd-hit.

        :param input_file: Ruta al archivo FASTA de entrada.
        :param output_prefix: Prefijo para los archivos de salida de cd-hit.
        :param threshold: Umbral de similitud para cd-hit. Por defecto es 0.7.
        """
        cd_hit_input = os.path.join(self.output_dir, f"{input_file}.fasta")
        cd_hit_output = os.path.join(self.output_dir, f"{input_file}")

        # Ejecutar cd-hit
        cd_hit(
            i=cd_hit_input,
            o=cd_hit_output,
            c=0.7,
            d=40,
            sc=1,
        )
        # # Leer y procesar el archivo de clústeres
        df_clstr = read_clstr(f"{cd_hit_output}.clstr")
        for _, row in df_clstr.iterrows():
            pdb_id = row["identifier"].split('_')[0]
            chain_number = row["identifier"].split('_')[1].split('|')[0]
            pdb_reference_id = self.session.query(PDBReference).filter_by(pdb_id=pdb_id).one().id
            chain_id = self.session.query(PDBChains).filter_by(pdb_reference_id=pdb_reference_id,
                                                               chain_number=chain_number).one().id
            cluster = Cluster(
                cluster_id=row['cluster'],
                pdb_chain_id=chain_id,
                is_representative=row['is_representative'],
                sequence_length=row['size'],
                identity=row['identity']
            )
            self.session.add(cluster)
        self.session.commit()
        logging.info(f"Clustering completado. Archivo de salida: {cd_hit_output}")


        return df_clstr

    def merge_fastas(self, pdb_ids, merge_name):
        """
        Combina archivos FASTA para un conjunto de IDs de PDB en un único
        archivo. Descarga los archivos FASTA que no están presentes en el
        directorio local.

        :param pdb_ids: Lista de identificadores de PDB para los cuales
            combinar los archivos FASTA.
        :type pdb_ids: list[str]
        """

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
