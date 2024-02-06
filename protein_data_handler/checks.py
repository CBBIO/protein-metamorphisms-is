from Bio.PDB import MMCIFParser, CEAligner

# Cargar las estructuras
parser = MMCIFParser()
structure_A = parser.get_structure("A", "./7A0C_B_0.cif")
structure_B = parser.get_structure("B", "./6R4Q_A_0.cif")

# Crear una instancia de CEAligner
aligner = CEAligner()

# Establecer la estructura de referencia
aligner.set_reference(structure_A)

# Alinear la estructura B con la estructura de referencia A
aligner.align(structure_A)

# Imprimir el RMSD
print(f"RMSD del alineamiento: {aligner.rms}")
