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
def test_cif_to_pdb_chain_id_modification(mock_pdbio, mock_mmcif_parser):
    # Crear mocks para las cadenas con diferentes IDs
    mock_chain_long_id = MagicMock(spec=['id'])
    mock_chain_long_id.id = "AB"  # Este ID debería ser cambiado a "X"
    mock_chain_short_id = MagicMock(spec=['id'])
    mock_chain_short_id.id = "A"  # Este ID debería permanecer igual

    # Simular la iteración sobre modelos y cadenas directamente
    mock_model = MagicMock()
    mock_model.__iter__.return_value = iter([mock_chain_long_id, mock_chain_short_id])

    mock_structure = MagicMock()
    mock_structure.__iter__.return_value = iter([mock_model])

    # Configurar el MMCIFParser mock para retornar la estructura mock
    mock_mmcif_parser.return_value.get_structure.return_value = mock_structure

    # Llamar a la función bajo prueba
    cif_to_pdb('dummy_path.cif', 'dummy_path.pdb')

    # Verificar las modificaciones de los IDs
    assert mock_chain_long_id.id == "X", "El ID de la cadena larga debe ser cambiado a 'X'"
    assert mock_chain_short_id.id == "A", "El ID de la cadena corta debe permanecer igual"

    # Verificar que set_structure y save se llamaron correctamente
    mock_pdbio().set_structure.assert_called_once_with(mock_structure)
    mock_pdbio().save.assert_called_once_with('dummy_path.pdb')