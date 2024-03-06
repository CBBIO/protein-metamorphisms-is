import multiprocessing
from datetime import timedelta, datetime
from unittest import TestCase, mock
from unittest.mock import patch, MagicMock, call, create_autospec
from multiprocessing.pool import ThreadPool as Pool
from sqlalchemy.orm.util import AliasedClass

from protein_metamorphisms_is.operations.structural_alignment import StructuralAlignmentManager
from protein_metamorphisms_is.sql.model import StructuralAlignmentQueue


class MockQueueItem:
    def __init__(self,queue_entry_id, alignment_type_id, align_task):
        self.queue_entry_id = queue_entry_id
        self.alignment_type_id = alignment_type_id
        self.align_task = align_task



def mock_align_task(item, conf):
    return f"Result for {item.queue_entry_id}"


class MockAlignTaskHolder:
    def __init__(self, align_task_func):
        self.align_task = align_task_func


class TestStructuralAlignment(TestCase):
    def setUp(self):
        # Configuration for CDHit instance
        self.conf = {
            "DB_USERNAME": "usuario",
            "DB_PASSWORD": "clave",
            "DB_HOST": "localhost",
            "DB_PORT": 5432,
            "DB_NAME": "BioData",
            'max_workers': 2,
            'constants': 'path/to/constants.yaml',
            'structural_alignment': {'batch_size': 1000, 'retry_count': 5, 'task_timeout': 50, 'types': ['1', '2']}
        }

        with patch('protein_metamorphisms_is.information_system.base.extractor.setup_logger'), \
                patch('protein_metamorphisms_is.sql.base.database_manager.create_engine'), \
                patch('protein_metamorphisms_is.sql.base.database_manager.sessionmaker', return_value=MagicMock()), \
                patch('protein_metamorphisms_is.operations.base.operator.handle_structural_alignment_types'), \
                patch('protein_metamorphisms_is.operations.base.operator.handle_structural_complexity_levels'), \
                patch('yaml.safe_load'), \
                patch('builtins.open', mock.mock_open(read_data="constants: value")):
            self.structural_alignment = StructuralAlignmentManager(self.conf)
            self.structural_alignment.types = {'1': MockAlignTaskHolder(mock_align_task)}

    @patch('protein_metamorphisms_is.operations.structural_alignment.StructuralAlignmentManager.execute_aligns')
    @patch('protein_metamorphisms_is.operations.structural_alignment.StructuralAlignmentManager.fetch_queue_items',
           return_value=[{'id': 1}])
    @patch('protein_metamorphisms_is.operations.structural_alignment.StructuralAlignmentManager.check_empty_queue',
           side_effect=[False, True])
    @patch('protein_metamorphisms_is.operations.structural_alignment.StructuralAlignmentManager.get_update_queue')
    @patch('protein_metamorphisms_is.operations.structural_alignment.StructuralAlignmentManager.fetch_tasks_info')
    def test_start(self, mock_fetch_tasks_info, mock_get_update_queue, mock_check_empty_queue, mock_fetch_queue_items,
                   mock_execute_aligns):
        self.structural_alignment.start()

        mock_fetch_tasks_info.assert_called_once()
        mock_get_update_queue.assert_called_once()
        mock_check_empty_queue.assert_called()
        mock_fetch_queue_items.assert_called_once()

    @patch('protein_metamorphisms_is.operations.structural_alignment.StructuralAlignmentManager.execute_aligns')
    @patch('protein_metamorphisms_is.operations.structural_alignment.StructuralAlignmentManager.fetch_queue_items',
           return_value=[{'id': 1}])
    @patch('protein_metamorphisms_is.operations.structural_alignment.StructuralAlignmentManager.check_empty_queue',
           side_effect=[False, True])
    @patch('protein_metamorphisms_is.operations.structural_alignment.StructuralAlignmentManager.get_update_queue')
    @patch('protein_metamorphisms_is.operations.structural_alignment.StructuralAlignmentManager.fetch_tasks_info')
    def test_start_with_exception(self, mock_fetch_tasks_info, mock_get_update_queue, mock_check_empty_queue,
                                  mock_fetch_queue_items, mock_execute_aligns):
        # Configura mock_fetch_tasks_info para lanzar una excepción
        mock_fetch_tasks_info.side_effect = Exception("Test exception")

        # Ejecuta el método bajo prueba y verifica que se maneje la excepción
        with self.assertRaises(Exception) as context:
            self.structural_alignment.start()

        # Verifica que el logger registró el error

        # Verifica que la excepción lanzada es la esperada

    @patch('importlib.import_module')
    def test_fetch_tasks_info(self, mock_import):
        # Configura el mock de la sesión para devolver tipos de alineación estructural simulados
        self.structural_alignment.session.query().all.return_value = [MagicMock(id='1', task_name='task1'),
                                                                      MagicMock(id='2', task_name='task2')]
        # Configura la configuración para incluir solo uno de los tipos de tareas

        self.structural_alignment.fetch_tasks_info()

        # Verifica que el módulo correcto se haya importado
        # Verifica que el tipo correcto se haya almacenado
        self.assertIn('1', self.structural_alignment.types)
        self.assertIn('2', self.structural_alignment.types)

    @patch('protein_metamorphisms_is.operations.structural_alignment.datetime', autospec=True)
    def test_get_update_queue(self, mock_datetime):
        # Configura la fecha/hora actual simulada
        now = datetime(2022, 1, 1, 12, 0, 0)
        mock_datetime.now.return_value = now

        # Simula las respuestas de las consultas a la base de datos
        self.structural_alignment.session.query().filter().all.side_effect = [
            # Primera llamada devuelve IDs de clusters ya en cola
            [('cluster1',), ('cluster2',)],
            # Segunda llamada devuelve clusters no en cola
            [MagicMock(id='cluster3'), MagicMock(id='cluster4')],
            # Tercera llamada devuelve entradas obsoletas
            [MagicMock(id=1, state=1, retry_count=0, updated_at=now - timedelta(days=1)),
             MagicMock(id=2, state=3, retry_count=1, updated_at=now - timedelta(days=1))]
        ]

        # Configura el diccionario de tipos para simular la existencia de un tipo de alineamiento

        self.structural_alignment.get_update_queue()

        # Verifica que se hayan añadido nuevas tareas a la cola
        # Aquí cambiamos el enfoque para verificar las propiedades de los objetos pasados a `session.add`
        added_calls = self.structural_alignment.session.add.call_args_list
        self.assertEqual(len(added_calls), 4)  # Asegúrate de que se haya llamado a `add` cuatro veces

        # Verifica las propiedades de las nuevas instancias de StructuralAlignmentQueue añadidas
        for call in added_calls[:2]:  # Asume que las primeras 2 llamadas son para añadir nuevas tareas
            queue_item = call[0][0]  # Accede al primer argumento de la llamada
            self.assertIsInstance(queue_item, StructuralAlignmentQueue)
            self.assertIn(queue_item.cluster_entry_id, ['cluster3', 'cluster4'])
            self.assertEqual(queue_item.state, 0)
            self.assertEqual(queue_item.alignment_type_id, '1')
            self.assertEqual(queue_item.retry_count, 0)

        # Para las entradas obsoletas, verifica que se hayan actualizado correctamente
        # Esto puede requerir acceder a las instancias de `StructuralAlignmentQueue` actualizadas y verificar sus propiedades

        # Verifica que se haya hecho commit a la sesión
        self.structural_alignment.session.commit.assert_called_once()

    def test_check_empty_queue(self):
        # Escenario 1: La cola está vacía
        self.structural_alignment.session.query().filter().count.return_value = 0
        self.assertTrue(self.structural_alignment.check_empty_queue())

        # Reinicia el mock para el siguiente escenario

        # Escenario 2: La cola tiene tareas pendientes o con errores
        self.structural_alignment.session.query().filter().count.return_value = 3
        self.assertFalse(self.structural_alignment.check_empty_queue())

    @patch('protein_metamorphisms_is.operations.structural_alignment.aliased')
    def test_fetch_queue_items(self, mock_aliased):
        # Crea un spec para el objeto que esperarías que sqlalchemy.orm.aliased devuelva
        # Esto puede ser un AliasedClass o cualquier clase que normalmente pasarías a aliased

        # Simula los resultados esperados de la consulta
        expected_queue_items = [
            MagicMock(id=1, alignment_type_id='1', cluster_id='cluster1', pdb_id='PDB1', chains='A', model=1,
                      rep_pdb_id='PDB1', rep_chains='A', rep_model=1),
            MagicMock(id=2, alignment_type_id='1', cluster_id='cluster2', pdb_id='PDB2', chains='B', model=1,
                      rep_pdb_id='PDB2', rep_chains='B', rep_model=1)
        ]
        self.structural_alignment.session.query().join().join().join().outerjoin().filter().order_by().limit().all.return_value = expected_queue_items

        # Llama al método que estás probando
        result = self.structural_alignment.fetch_queue_items()

        # Verifica que el resultado sea como se espera
        self.assertEqual(result, expected_queue_items)
        # Asegúrate de que la consulta se construyó correctamente
        self.structural_alignment.session.query().join().join().join().outerjoin().filter().order_by().limit.assert_called_with(
            self.conf['structural_alignment']['batch_size'])

    @patch('protein_metamorphisms_is.operations.structural_alignment.StructuralAlignmentManager.insert_results')
    @patch('multiprocessing.pool.ThreadPool', autospec=True)
    def test_execute_aligns(self, mock_insert_results, mock_pool):
        # Simula los elementos de la cola
        queue_items = [MockQueueItem(1, '1', mock_align_task), MockQueueItem(2, '1', mock_align_task)]

        # Ejecutar el método bajo prueba
        self.structural_alignment.execute_aligns(queue_items)

        # Verificaciones
        self.structural_alignment.insert_results.assert_called_once()

    @patch('multiprocessing.pool.ApplyResult', autospec=True)
    def test_execute_aligns_with_timeout_error(self, mock_pool):
        # Configurar el mock de ThreadPool para simular un TimeoutError en apply_async.get
        self.structural_alignment.insert_results = MagicMock()
        mock_task_get = MagicMock(side_effect=multiprocessing.TimeoutError)
        mock_task = MagicMock()
        mock_task.get = mock_task_get
        mock_pool.return_value = mock_task

        # Simular elementos de la cola
        queue_items = [MockQueueItem(1, '1', mock_align_task), MockQueueItem(2, '1', mock_align_task)]

        # Ejecutar el método bajo prueba
        self.structural_alignment.execute_aligns(queue_items)



    @patch('protein_metamorphisms_is.operations.structural_alignment.StructuralAlignmentResults')
    def test_insert_results(self, mock_StructuralAlignmentResults):
        # Simular resultados de tareas de alineación estructural
        successful_result = (1, {'alignment_score': 100})
        failed_result = (2, {'error_message': 'Error processing alignment'})

        # Crear instancias mock de StructuralAlignmentQueue para simular las entradas de la cola
        mock_queue_item_success = MagicMock(spec=StructuralAlignmentQueue)
        mock_queue_item_failed = MagicMock(spec=StructuralAlignmentQueue)
        self.structural_alignment.session.query.return_value.filter_by.return_value.first.side_effect = [
            mock_queue_item_success, mock_queue_item_failed]

        # Ejecutar el método bajo prueba
        self.structural_alignment.insert_results([successful_result, failed_result])

        # Verificaciones para el resultado exitoso
        mock_StructuralAlignmentResults.assert_called_once_with(alignment_score=100)
        self.structural_alignment.session.add.assert_any_call(mock_StructuralAlignmentResults.return_value)
        self.assertEqual(mock_queue_item_success.state, 2)
        self.assertEqual(mock_queue_item_failed.state, 3)
        self.assertEqual(mock_queue_item_failed.error_message, 'Error processing alignment')
        self.structural_alignment.session.commit.assert_called_once()
