from Bio import SeqIO
import pandas as pd

# Ruta al archivo FASTA
fasta_file = 'data/uniprot_sprot.fasta'

# Lista para almacenar los IDs
ids = []

# Leer el archivo FASTA usando BioPython
for record in SeqIO.parse(fasta_file, 'fasta'):
    ids.append(record.id.split('|')[1])

# Crear un DataFrame de pandas sin índice
df = pd.DataFrame(ids, columns=['ID'])

# Guardar el DataFrame en un archivo CSV sin índice
csv_file = 'data/goa_ids.csv'
df.to_csv(csv_file, index=False)


print(f"IDs extraídos y guardados en {csv_file}")
