import re
from Bio import SeqIO
from Bio.PDB import MMCIFParser, PDBIO
import gemmi


def extract_float(s):
    """
    Extrae el primer número flotante de una cadena de texto.

    Esta función busca en la cadena proporcionada el primer patrón que
    corresponda a un número flotante. Si lo encuentra, lo devuelve como
    un valor de tipo float. Si no encuentra ningún número flotante,
    devuelve None.

    :param s: La cadena de texto de la cual extraer el número flotante.
    :type s: str
    :return: El primer número flotante encontrado en la cadena, o None si no
             se encuentra ninguno.
    :rtype: float or None

    :Example:

    >>> extract_float("1.75 A")
    1.75
    >>> extract_float("El valor es -3.5")
    -3.5
    >>> extract_float("-")
    None

    .. note::
        Esta función solo extrae el primer número flotante encontrado
        en la cadena.
    """

    match = re.search(r"(\d+\.\d+)", s)
    if match:
        return float(match.group(1))
    else:
        return None


def safe_convert_to_int(value):
    """
    Intenta convertir un valor a un entero.

    Esta función intenta convertir el valor dado a un entero. Si la conversión
    falla debido a un ValueError (por ejemplo, si el valor es una cadena vacía
    o una cadena que no representa un entero), maneja la excepción retornando
    un valor por defecto.

    :param value: El valor que se intentará convertir a entero.
        Puede ser de cualquier tipo que sea convertible a entero,
        típicamente una cadena o un número.
    :type value: str o int
    :return: El valor convertido a entero, o None si la conversión falla.
    :rtype: int o None

    """
    try:
        return int(value)
    except ValueError:
        # Manejar el error o retornar un valor por defecto como None o 0
        # dependiendo de tus requerimientos
        return None


def process_chain_string(chain):
    """
    Procesa una cadena que representa un objeto de cadena y su rango.

    Esta función toma una cadena en el formato
    'nombreCadena=rangoInicio-rangoFin',
    la divide y convierte los rangos de inicio y fin a enteros.
    Retorna el nombre de la cadena y los valores de inicio y fin
    del rango como enteros.

    :param chain_obj: La cadena de objeto de cadena a procesar.
    :type chain_obj: str
    :return: El nombre de la cadena, el valor de inicio y
    el valor de fin del rango.

    """

    chain = chain.replace(" ", "")
    if '=' not in chain:
        return None, None, None

    chain_name, seq_range = chain.split('=')
    try:
        start, end = seq_range.split('-')
        start = safe_convert_to_int(start)
        end = safe_convert_to_int(end)
    except Exception:
        start, end = None, None

    return chain_name, start, end


def extract_and_parse_fasta(file_path):
    """
    Extrae y analiza las secuencias de un archivo FASTA.

    :param file_path: Ruta del archivo FASTA.
    :return: Lista de tuplas con (pdb_id, chain, sequence).
    """
    sequences = []

    with open(file_path, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            header = record.description
            # Extraer los campos del encabezado
            parts = header.split('|')
            pdb_id = parts[0].split('_')[0]
            chain_number = parts[0].split('_')[1]

            chain = parts[1].replace('Chains ', '').replace('Chain ', '').replace(',', '/').replace(' ', '')
            chain = auth_chain_mapping(chain)
            sequence = str(record.seq)
            # Agregar la tupla a la lista
            sequences.append((pdb_id, chain_number, chain, sequence))

    return sequences


def auth_chain_mapping(chain):
    elementos = chain.split('/')
    elementos_transformados = []
    for elemento in elementos:
        match = re.search(r'[A-Za-z]+\[auth([A-Za-z0-9]+)\]', elemento)
        if match:
            elementos_transformados.append(match.group(1))
        else:
            elementos_transformados.append(elemento)

    resultado = '/'.join(elementos_transformados)
    return resultado


def cif_to_pdb(cif_path, pdb_path):
    """Convierte un archivo CIF a PDB."""

    parser = MMCIFParser()
    structure = parser.get_structure('ID', cif_path)
    io = PDBIO()
    io.set_structure(structure)

    for model in structure:
        for chain in model:
            if len(chain.id) > 1:
                chain.id = "X"
    io.save(pdb_path)


def get_chain_to_accession_map(cif_path):
    """
    Extracts the mapping of chain names to UniProt accessions from a CIF file.

    Args:
        cif_path (str): Path to the CIF file.

    Returns:
        dict: A dictionary mapping chain names to UniProt accessions.
    """
    doc = gemmi.cif.read(cif_path)
    block = doc.sole_block()

    # Initialize mappings
    chain_accession_map = {}

    try:
        # Extract _entity_poly fields
        entity_ids_poly = block.find_values('_entity_poly.entity_id')
        pdbx_strand_ids = block.find_values('_entity_poly.pdbx_strand_id')

        # Extract _struct_ref fields
        db_names = block.find_values('_struct_ref.db_name')
        entity_ids_ref = block.find_values('_struct_ref.entity_id')
        pdbx_db_accessions = block.find_values('_struct_ref.pdbx_db_accession')

        # Map entity IDs to UniProt accessions
        struct_ref_map = {
            entity_id: accession
            for db_name, entity_id, accession in zip(db_names, entity_ids_ref, pdbx_db_accessions)
            if db_name == 'UNP' and accession  # Only consider UniProt entries
        }

        # Map entity IDs to chain names
        entity_chain_map = {
            entity_id: chains.split(',')
            for entity_id, chains in zip(entity_ids_poly, pdbx_strand_ids)
            if chains
        }

        # Combine mappings to produce chain -> accession map
        for entity_id, chains in entity_chain_map.items():
            if entity_id in struct_ref_map:
                accession = struct_ref_map[entity_id]
                for chain in chains:
                    chain_accession_map[chain.strip()] = accession
    except KeyError as e:
        print(f"Missing key in CIF file: {e}")
    return chain_accession_map
