import unittest
from protein_data_handler.helpers.parser.parser import extract_float, procesar_chain_string


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


if __name__ == '__main__':
    unittest.main()
