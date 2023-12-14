import unittest

from Bio import SeqRecord, Seq

from protein_data_handler.helpers.parser.parser import extract_float, procesar_chain_string, extract_and_parse_fasta, \
    auth_chain_mapping
from unittest.mock import patch, mock_open

class TestExtractFloat(unittest.TestCase):
    def test_extract_float_with_valid_float(self):
        self.assertEqual(extract_float("1.75 A"), 1.75)
        self.assertEqual(extract_float("0.99"), 0.99)
        self.assertEqual(extract_float("123.456"), 123.456)

    def test_extract_float_with_no_float(self):
        self.assertIsNone(extract_float("-"))
        self.assertIsNone(extract_float("ABC"))
        self.assertIsNone(extract_float("123"))

    def test_cadena_valida(self):
        resultado = procesar_chain_string("A=100-200")
        self.assertEqual(resultado, ("A", 100, 200))

    def test_cadena_sin_guion(self):
        resultado = procesar_chain_string("A=100200")
        self.assertEqual(resultado, ("A", None, None))

    def test_cadena_sin_igual(self):
        resultado = procesar_chain_string("A100-200")
        self.assertEqual(resultado, (None, None, None))

    def test_cadena_con_valores_no_numericos(self):
        resultado = procesar_chain_string("A=abc-def")
        self.assertEqual(resultado, ("A", None, None))

    def test_cadena_vacia(self):
        resultado = procesar_chain_string("")
        self.assertEqual(resultado, (None, None, None))

    def test_cadena_con_espacios(self):
        resultado = procesar_chain_string(" A = 100 - 200 ")
        self.assertEqual(resultado, ("A", 100, 200))

    def test_cadena_con_espacios_y_formato_incorrecto(self):
        resultado = procesar_chain_string(" A = 100200 ")
        self.assertEqual(resultado, ("A", None, None))

    @patch('protein_data_handler.helpers.parser.parser.open', new_callable=mock_open, read_data=">PDB1_A\nATCG")
    @patch('protein_data_handler.helpers.parser.parser.SeqIO.parse')
    def test_extract_and_parse_fasta(self, mock_seqio_parse, mock_file):
        mock_seqio_parse.return_value = [SeqRecord.SeqRecord(Seq.Seq("ATCG"), id="PDB1_A|Chain A", description="PDB1_A|Chain A")]

        file_path = 'fake_path.fasta'
        sequences = extract_and_parse_fasta(file_path)

        self.assertEqual(sequences, [("PDB1", "A", "A", "ATCG")])

        mock_file.assert_called_with(file_path, 'r')
        mock_seqio_parse.assert_called()

    def test_auth_chain_mapping(self):
        self.assertEqual(auth_chain_mapping("ChainA[authA]/ChainB[authB]"), "A/B")
        self.assertEqual(auth_chain_mapping("ChainA/ChainB"), "ChainA/ChainB")
        self.assertEqual(auth_chain_mapping("Chain[authA]"), "A")
        self.assertEqual(auth_chain_mapping(""), "")
if __name__ == '__main__':
    unittest.main()
