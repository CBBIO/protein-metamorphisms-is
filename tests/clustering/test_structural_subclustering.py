import os
import unittest

import pytest
from sqlalchemy.dialects.mssql.information_schema import sequences

from protein_metamorphisms_is.helpers.config.yaml import read_yaml_config
from protein_metamorphisms_is.operation.clustering.sequence_clustering import SequenceClustering
from protein_metamorphisms_is.operation.clustering.structural_subclustering import StructuralSubClustering
from protein_metamorphisms_is.operation.embedding.sequence_embedding import SequenceEmbeddingManager
from protein_metamorphisms_is.operation.extraction.uniprot import UniProtExtractor
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


@pytest.mark.order(8)
class TestStructuralSubClustering(unittest.TestCase):

    def setUp(self):
        config_path = self.get_config_path()
        config = read_yaml_config(config_path)
        self.clustering = StructuralSubClustering(config)

    def get_config_path(self):
        """Devuelve la ruta al archivo de configuración."""
        return os.path.join('tests/config/', "config.yaml")  # Ruta del archivo YAML de prueba

    def test_entities_created(self):
        """Verifica que se crean instancias en todas las entidades tras ejecutar el subclustering."""
        # Ejecutar el proceso de subclustering
        self.clustering.start()

        # Verificar que los clusters están presentes en la base de datos
        clusters_count = self.clustering.session.query(Cluster).count()
        self.assertGreater(clusters_count, 0, "No se generaron clusters en la base de datos.")

        # Verificar que los subclusters están presentes en la base de datos
        subclusters_count = self.clustering.session.query(Subcluster).count()
        self.assertGreater(subclusters_count, 0, "No se generaron subclusters en la base de datos.")

        # Verificar que las entradas de los subclusters están presentes
        subcluster_entries_count = self.clustering.session.query(SubclusterEntry).count()
        self.assertGreater(subcluster_entries_count, 0, "No se generaron entradas de subclusters.")

        # Verificar que cada subcluster tiene al menos una entrada asociada
        subclusters = self.clustering.session.query(Subcluster).all()
        for subcluster in subclusters:
            subcluster_entries = self.clustering.session.query(SubclusterEntry).filter_by(
                subcluster_id=subcluster.id).count()
            self.assertGreater(subcluster_entries, 0, f"El subcluster {subcluster.id} no tiene entradas asociadas.")

        # Verificar que todas las entradas de los subclusters tienen un embedding válido asociado
        subcluster_entries = self.clustering.session.query(SubclusterEntry).all()
        for entry in subcluster_entries:
            embedding = self.clustering.session.query(Structure3Di).get(entry.structure_3di_id)
            self.assertIsNotNone(embedding, f"El embedding con ID {entry.structure_3di_id} no existe.")

    if __name__ == '__main__':
        unittest.main()