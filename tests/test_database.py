import unittest
from unittest.mock import patch, MagicMock

from protein_data_handler.helpers.database.database import create_session


class TestCreateSession(unittest.TestCase):
    @patch("protein_data_handler.helpers.database.database.create_engine")
    @patch("protein_data_handler.helpers.database.database.sessionmaker")
    def test_creacion_correcta_sesion(
        self, mock_sessionmaker, mock_create_engine
    ):
        dummy_uri = "sqlite:///test.db"
        mock_session = MagicMock()
        mock_sessionmaker.return_value = MagicMock(return_value=mock_session)

        session = create_session(dummy_uri)

        mock_create_engine.assert_called_with(dummy_uri)
        mock_sessionmaker.assert_called_with(
            bind=mock_create_engine.return_value
        )
        self.assertEqual(session, mock_session)


if __name__ == "__main__":
    unittest.main()
