import unittest
from protein_data_handler.helpers.parser.parser import extract_float

class TestExtractFloat(unittest.TestCase):
    def test_extract_float_with_valid_float(self):
        self.assertEqual(extract_float("1.75 A"), 1.75)
        self.assertEqual(extract_float("0.99"), 0.99)
        self.assertEqual(extract_float("123.456"), 123.456)

    def test_extract_float_with_no_float(self):
        self.assertIsNone(extract_float("-"))
        self.assertIsNone(extract_float("ABC"))
        self.assertIsNone(extract_float("123"))


if __name__ == '__main__':
    unittest.main()
