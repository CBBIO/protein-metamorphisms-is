import os
import tempfile
from concurrent.futures import ThreadPoolExecutor
from io import StringIO
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MafftCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from sqlalchemy import func

from protein_data_handler.models.uniprot import PDBChain, UniprotChain, PDBReference, UniProtPDBAlignment

import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


class UniProtPDBMapping:
    """
    Clase para mapear y alinear secuencias de UniProt y PDB.

    :param session: Sesión de SQLAlchemy para interactuar con la base de datos.
    """

    def __init__(self, session):
        self.session = session

    def realizar_consulta_cadenas_iguales(self):
        """
        Realiza una consulta a la base de datos para encontrar cadenas iguales en UniProt y PDB.

        :return: Resultado de la consulta con secuencias y metadatos.
        """
        logging.info("Realizando consulta para encontrar cadenas iguales en UniProt y PDB.")
        try:
            result = self.session.query(
                UniprotChain.sequence.label('uniprot_sequence'),
                PDBChain.sequence.label('pdb_sequence'),
                UniprotChain.chain.label('uniprot_chain'),
                PDBReference.pdb_id,
                PDBChain.chain.label('pdb_chain'),
                func.length(PDBChain.sequence).label('length_pdb_sequence'),
                func.length(UniprotChain.sequence).label('length_uniprot_sequence'),
                PDBReference.id.label('pdb_reference_id')
            ).join(
                PDBReference, PDBChain.pdb_reference_id == PDBReference.id
            ).join(
                UniprotChain, PDBReference.id == UniprotChain.pdb_reference_id
            ).filter(
                PDBChain.chain == UniprotChain.chain
            ).all()
            logging.info(f"Consulta completada con éxito. Número de registros encontrados: {len(result)}")
            return result
        finally:
            self.session.close()
            logging.info("Sesión de base de datos cerrada.")

    def volcar_datos_alineamiento(self, pares):
        """
        Procesa y almacena los datos de alineamiento en la base de datos utilizando múltiples hilos.

        :param pares: Lista de pares de secuencias para alinear y almacenar.
        """
        logging.info("Iniciando el proceso de volcado de datos de alineamiento en múltiples hilos.")
        with ThreadPoolExecutor() as executor:
            executor.map(self.procesar_par, pares)

        self.session.commit()
        logging.info("Datos de alineamiento almacenados con éxito en la base de datos.")

    def procesar_par(self, par):
        """
        Procesa un par de secuencias para alinear y almacenar en la base de datos.

        :param par: Tupla con las secuencias y metadatos a procesar.
        """
        logging.info(f"Procesando par: {par}")
        uniprot_seq, pdb_seq, chain, pdb_id = par[:4]
        pdb_reference_id = self.session.query(PDBReference.id).filter_by(pdb_id=pdb_id).scalar()
        uniprot_sequence_aligned, pdb_sequence_aligned = self.alinear_secuencias_mafft(uniprot_seq, pdb_seq)

        alignment = UniProtPDBAlignment(
            chain=chain,
            pdb_reference_id=pdb_reference_id,
            uniprot_sequence_aligned=uniprot_sequence_aligned,
            pdb_sequence_aligned=pdb_sequence_aligned
        )
        self.session.add(alignment)

    def alinear_secuencias_mafft(self, seq1, seq2):
        """
        Alinea dos secuencias utilizando MAFFT.

        :param seq1: Primera secuencia para alinear.
        :param seq2: Segunda secuencia para alinear.
        :return: Tupla con las secuencias alineadas.
        """
        logging.info(f"Iniciando alineación MAFFT para las secuencias: {seq1[:10]}..., {seq2[:10]}...")

        # Crear un archivo temporal para las secuencias
        with tempfile.NamedTemporaryFile(mode='w+', suffix='.fasta', delete=False) as temp_file:
            record1 = SeqRecord(Seq(seq1), id="Seq1")
            record2 = SeqRecord(Seq(seq2), id="Seq2")
            SeqIO.write([record1, record2], temp_file, "fasta")
            temp_file_path = temp_file.name

        # Ejecutar MAFFT utilizando el archivo temporal
        mafft_cline = MafftCommandline(input=temp_file_path)
        stdout, stderr = mafft_cline()

        # Leer el resultado del alineamiento
        align = AlignIO.read(StringIO(stdout), "fasta")
        logging.info("Alineación completada con éxito.")

        # Eliminar el archivo temporal
        os.remove(temp_file_path)

        return str(align[0].seq), str(align[1].seq)
