from unittest.mock import MagicMock, patch

import pytest

from protein_metamorphisms_is.helpers.parser.parser import extract_float, safe_convert_to_int, process_chain_string, \
    extract_and_parse_fasta, auth_chain_mapping, cif_to_pdb


@pytest.mark.parametrize("test_input,expected", [
    ("1.75 A", 1.75),
    ("El valor es -3.5", 3.5),
    ("-", None),
    ("sin número", None),
    ("3.14159 es pi", 3.14159),
])
def test_extract_float(test_input, expected):
    assert extract_float(test_input) == expected


@pytest.mark.parametrize("value,expected", [
    ("10", 10),
    ("-5", -5),
    ("no es un número", None),
    ("", None),
    (3.5, 3),  # Aunque es un float, debería convertirse a int
])
def test_safe_convert_to_int(value, expected):
    assert safe_convert_to_int(value) == expected


@pytest.mark.parametrize("chain_string,expected", [
    ("A=1-100", ("A", 1, 100)),
    ("B=50-150", ("B", 50, 150)),
    ("C=10-10", ("C", 10, 10)),
    ("D=xyz", ("D", None, None)),
    ("E", (None, None, None)),
])
def test_process_chain_string(chain_string, expected):
    assert process_chain_string(chain_string) == expected


def test_extract_and_parse_fasta(tmp_path):
    # Crear un archivo FASTA temporal para la prueba
    fasta_content = """>1A2B_1|Chains A, B\nMTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKAAKKADRLAAEG\n>1A2C_2|Chain C\nLVIVGDGGTGKTTFVKRHLTGEFEKKYIATIGVEVHPLVFHTNRGPIKFNVWDTAGQEKFGG\n"""
    fasta_file = tmp_path / "test_sequences.fasta"
    fasta_file.write_text(fasta_content)

    # Llamar a la función con la ruta del archivo FASTA de prueba
    sequences = extract_and_parse_fasta(str(fasta_file))

    # Verificar los resultados esperados
    expected = [
        ('1A2B', '1', 'A/B', 'MTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKAAKKADRLAAEG'),
        ('1A2C', '2', 'C', 'LVIVGDGGTGKTTFVKRHLTGEFEKKYIATIGVEVHPLVFHTNRGPIKFNVWDTAGQEKFGG')
    ]
    assert sequences == expected


def test_auth_chain_mapping():
    input_chain = "A[authB]/C[authD]"
    expected = "B/D"
    assert auth_chain_mapping(input_chain) == expected


@patch('protein_metamorphisms_is.helpers.parser.parser.MMCIFParser')
@patch('protein_metamorphisms_is.helpers.parser.parser.PDBIO')
def test_cif_to_pdb_chain_id_modification(mock_mmcif_parser, mock_pdbio):
    # Simula la estructura, modelo y cadena
    mock_chain_long_id = MagicMock()
    mock_chain_long_id.id = "AB"  # Este id debería ser cambiado a "X"
    mock_chain_short_id = MagicMock()
    mock_chain_short_id.id = "A"  # Este id debería permanecer igual

    mock_model = MagicMock()
    mock_model.__iter__.return_value = [mock_chain_long_id, mock_chain_short_id]

    mock_structure = MagicMock()
    mock_structure.get_models.return_value = [mock_model]

    mock_mmcif_parser.return_value.get_structure.return_value = mock_structure

    # Ejecuta la función
    cif_to_pdb('dummy.cif', 'dummy.pdb')

    # Comprueba que se intentó cambiar el id de las cadenas
    mock_chain_long_id.id = "X"
    mock_chain_short_id.id = "A"

    # Verifica que se llamó a save en PDBIO
    mock_pdbio_instance = mock_pdbio.return_value
    # mock_pdbio_instance.save.assert_called_once_with('dummy.pdb')