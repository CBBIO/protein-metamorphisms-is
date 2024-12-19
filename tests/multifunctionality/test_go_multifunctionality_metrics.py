import os
import unittest

import pytest
from sqlalchemy.dialects.mssql.information_schema import sequences

from protein_metamorphisms_is.helpers.config.yaml import read_yaml_config
from protein_metamorphisms_is.operation.clustering.sequence_clustering import SequenceClustering
from protein_metamorphisms_is.operation.clustering.structural_subclustering import StructuralSubClustering
from protein_metamorphisms_is.operation.embedding.sequence_embedding import SequenceEmbeddingManager
from protein_metamorphisms_is.operation.extraction.uniprot import UniProtExtractor


from protein_metamorphisms_is.operation.functional.multifunctionality.go_multifunctionality_metrics import \
    GoMultifunctionalityMetrics
from protein_metamorphisms_is.operation.structural_alignment.structural_alignment import StructuralAlignmentManager
from protein_metamorphisms_is.sql.model.entities.embedding.sequence_embedding import SequenceEmbedding, \
    SequenceEmbeddingType
from protein_metamorphisms_is.sql.model.entities.embedding.structure_3di import Structure3Di
from protein_metamorphisms_is.sql.model.entities.go_annotation.go_annotation import ProteinGOTermAnnotation
from protein_metamorphisms_is.sql.model.entities.go_annotation.go_term import GOTerm
from protein_metamorphisms_is.sql.model.entities.protein.accesion import Accession
from protein_metamorphisms_is.sql.model.entities.protein.protein import Protein
from protein_metamorphisms_is.sql.model.entities.sequence.sequence import Sequence
from protein_metamorphisms_is.sql.model.entities.structure.structure import Structure
from protein_metamorphisms_is.sql.model.operational.clustering.cluster import Cluster, Subcluster, SubclusterEntry


@pytest.mark.order(10)
class TestGoMultifunctionalityMetrics(unittest.TestCase):

    def setUp(self):
        config_path = self.get_config_path()
        config = read_yaml_config(config_path)
        self.multifunctionality = GoMultifunctionalityMetrics(config)

    def get_config_path(self):
        """Devuelve la ruta al archivo de configuración."""
        return os.path.join('tests/config/', "config.yaml")  # Ruta del archivo YAML de prueba

    def test_entities_created(self):
        """Verifica que se crean instancias en todas las entidades tras ejecutar el subclustering."""
        # Ejecutar el proceso de subclusteringhe a
        self.multifunctionality.start()



    if __name__ == '__main__':
        unittest.main()