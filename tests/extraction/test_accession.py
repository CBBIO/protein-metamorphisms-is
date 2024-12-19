import re
import unittest
import os

import pytest

from protein_metamorphisms_is.helpers.config.yaml import read_yaml_config
from protein_metamorphisms_is.operation.extraction.accessions import AccessionManager
from protein_metamorphisms_is.operation.extraction.pdb import PDBExtractor
from protein_metamorphisms_is.operation.extraction.uniprot import UniProtExtractor
from protein_metamorphisms_is.operation.embedding.sequence_embedding import SequenceEmbeddingManager
from protein_metamorphisms_is.operation.embedding.structure_3di import Structure3DiManager
from protein_metamorphisms_is.operation.clustering.sequence_clustering import SequenceClustering
from protein_metamorphisms_is.operation.clustering.structural_subclustering import StructuralSubClustering
from protein_metamorphisms_is.operation.functional.annotation_transfer.sequence_go_annotation import \
    SequenceGOAnnotation
from protein_metamorphisms_is.operation.structural_alignment.structural_alignment import StructuralAlignmentManager
from protein_metamorphisms_is.operation.functional.multifunctionality.go_multifunctionality_metrics import \
    GoMultifunctionalityMetrics
from protein_metamorphisms_is.sql.model.entities.protein.accesion import Accession

@pytest.mark.order(1)
class TestAccessionManager(unittest.TestCase):

    def setUp(self):
        config_path = self.get_config_path()
        config = read_yaml_config(config_path)
        self.extractor = AccessionManager(config)

    def get_config_path(self):
        """Devuelve la ruta al archivo de configuración."""
        return os.path.join('tests/config/', "config.yaml")  # Ruta del archivo YAML de prueba

    def test_conf(self):
        """Verifica que la configuración se haya cargado correctamente."""
        self.assertIn("DB_USERNAME", self.extractor.conf)  # Ajusta según tu YAML
        self.assertIn("DB_HOST", self.extractor.conf)  # Ajusta según tu YAML
        self.assertIn("DB_NAME", self.extractor.conf)  # Ajusta según tu YAML

    def test_api_results(self):
        """Verifica que se carguen accesiones desde la API."""
        self.extractor.fetch_accessions_from_api()

        # Verificamos que se hayan cargado accesiones
        accessions = self.extractor.session.query(Accession).all()
        self.assertGreater(len(accessions), 0, "No se cargaron accesiones desde la API")

    def test_csv_results(self):
        """Verifica que se carguen accesiones desde la API."""
        self.extractor.load_accessions_from_csv()

    def test_unique_accessions(self):
        """Verifica que las accesiones cargadas sean únicas."""
        accessions = self.extractor.session.query(Accession).all()
        accession_codes = [acc.code for acc in accessions]

        # Verificamos que no haya accesiones duplicadas
        self.assertEqual(len(accession_codes), len(set(accession_codes)), "Hay accesiones duplicadas en los resultados")


    def test_valid_accession_format(self):
        """Verifica que las accesiones tengan un formato válido (de 6 a 10 caracteres alfanuméricos)."""
        self.extractor.fetch_accessions_from_api()

        accessions = self.extractor.session.query(Accession).all()

        # Definimos un patrón para los códigos de acceso: 6 a 10 caracteres alfanuméricos
        accession_pattern = re.compile(r'^[A-Za-z0-9]{6,10}$')

        for accession in accessions:
            self.assertTrue(accession_pattern.match(accession.code),
                            f"Formato inválido para la accesión {accession.code}")


