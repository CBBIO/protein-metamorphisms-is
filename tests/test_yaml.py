import unittest
from unittest.mock import patch, mock_open

from protein_metamorphisms_is.helpers.config.yaml import read_yaml_config


class TestReadYamlConfig(unittest.TestCase):
    @patch("builtins.open", new_callable=mock_open, read_data="key: value")
    @patch("yaml.safe_load")
    def test_lectura_correcta_archivo_yaml(self, mock_safe_load, mock_file):
        mock_safe_load.return_value = {"key": "value"}

        result = read_yaml_config("dummy_file_path.yaml")

        mock_file.assert_called_with("dummy_file_path.yaml", "r")
        mock_safe_load.assert_called()
        self.assertEqual(result, {"key": "value"})

    @patch("builtins.open", new_callable=mock_open, read_data="key: value")
    def test_manejo_excepciones_archivo_yaml(self, mock_file):
        mock_file.side_effect = FileNotFoundError

        with self.assertRaises(FileNotFoundError):
            read_yaml_config("nonexistent_file_path.yaml")


if __name__ == "__main__":
    unittest.main()
