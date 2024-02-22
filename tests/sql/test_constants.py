import pytest
from unittest.mock import MagicMock
from protein_metamorphisms_is.sql.model import StructuralComplexityLevel, StructuralAlignmentType
from protein_metamorphisms_is.sql.constants import handle_structural_complexity_levels, handle_structural_alignment_types

def test_handle_structural_complexity_levels():
    # Simular la sesión y los datos de entrada
    mock_session = MagicMock()
    mock_query = mock_session.query.return_value
    mock_filter_by = mock_query.filter_by.return_value
    constants = {
        'structural_complexity_levels': [
            {'name': 'Low', 'description': 'Low complexity'},
            {'name': 'High', 'description': 'High complexity'}
        ]
    }

    # Configurar el mock para simular que no existen entradas previas en la base de datos
    mock_filter_by.first.side_effect = [None, None]

    # Llamar a la función bajo prueba
    handle_structural_complexity_levels(mock_session, constants)

    # Verificar que se hayan creado y añadido los nuevos niveles de complejidad
    assert mock_session.add.call_count == len(constants['structural_complexity_levels'])
    mock_session.commit.assert_called_once()

def test_handle_structural_alignment_types():
    # Simular la sesión y los datos de entrada
    mock_session = MagicMock()
    mock_query = mock_session.query.return_value
    mock_filter_by = mock_query.filter_by.return_value
    constants = {
        'structural_alignment_types': [
            {'name': 'FATCAT', 'description': 'Flexible structure alignment'},
            {'name': 'CE', 'description': 'Combinatorial Extension'}
        ]
    }

    # Configurar el mock para simular que no existen entradas previas en la base de datos
    mock_filter_by.first.side_effect = [None, None]

    # Llamar a la función bajo prueba
    handle_structural_alignment_types(mock_session, constants)

    # Verificar que se hayan creado y añadido los nuevos tipos de alineamiento
    assert mock_session.add.call_count == len(constants['structural_alignment_types'])
    mock_session.commit.assert_called_once()
