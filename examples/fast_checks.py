from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MafftCommandline
from io import StringIO

def align_sequences_with_mafft(seq1, seq2):
    # Crear registros de secuencia
    record1 = SeqRecord(Seq(seq1), id="Seq1")
    record2 = SeqRecord(Seq(seq2), id="Seq2")

    # Guardar en un archivo FASTA
    SeqIO.write([record1, record2], "temp_sequences.fasta", "fasta")

    # Crear la línea de comando para MAFFT
    mafft_cline = MafftCommandline(input="temp_sequences.fasta")

    # Ejecutar MAFFT y capturar la salida
    stdout, stderr = mafft_cline()

    # Leer el alineamiento desde la salida estándar
    align = AlignIO.read(StringIO(stdout), "fasta")

    # Convertir las secuencias alineadas a strings
    aligned_seq1 = str(align[0].seq)
    aligned_seq2 = str(align[1].seq)

    return aligned_seq1, aligned_seq2

# Ejemplo de uso
seq1 = "ATCGTAC"
seq2 = "ATGAC"
aligned_seq1, aligned_seq2 = align_sequences_with_mafft(seq1, seq2)
print("Secuencia 1 Alineada:", aligned_seq1)
print("Secuencia 2 Alineada:", aligned_seq2)
