import os
import unittest

import pytest
from sqlalchemy.dialects.mssql.information_schema import sequences

from protein_metamorphisms_is.helpers.config.yaml import read_yaml_config
from protein_metamorphisms_is.operation.clustering.sequence_clustering import SequenceClustering
from protein_metamorphisms_is.operation.embedding.sequence_embedding import SequenceEmbeddingManager
from protein_metamorphisms_is.operation.extraction.uniprot import UniProtExtractor
from protein_metamorphisms_is.sql.model.entities.embedding.sequence_embedding import SequenceEmbedding, \
    SequenceEmbeddingType
from protein_metamorphisms_is.sql.model.entities.go_annotation.go_annotation import ProteinGOTermAnnotation
from protein_metamorphisms_is.sql.model.entities.go_annotation.go_term import GOTerm
from protein_metamorphisms_is.sql.model.entities.protein.accesion import Accession
from protein_metamorphisms_is.sql.model.entities.protein.protein import Protein
from protein_metamorphisms_is.sql.model.entities.sequence.sequence import Sequence
from protein_metamorphisms_is.sql.model.entities.structure.structure import Structure
from protein_metamorphisms_is.sql.model.operational.clustering.cluster import ClusterEntry, Cluster


@pytest.mark.order(7)
class TestSequenceClustering(unittest.TestCase):

    def setUp(self):
        config_path = self.get_config_path()
        config = read_yaml_config(config_path)
        self.clustering = SequenceClustering(config)

    def get_config_path(self):
        """Devuelve la ruta al archivo de configuración."""
        return os.path.join('tests/config/', "config.yaml")  # Ruta del archivo YAML de prueba

    def test_entities_created(self):
        """Verifica que se crean instancias en todas las entidades tras ejecutar el embedder."""
        # Ejecutar el manager para procesar embeddings
        self.clustering.start()

        # Verificar que las secuencias están presentes en la base de datos
        sequences_count = self.clustering.session.query(Sequence).count()
        self.assertGreater(sequences_count, 0, "No se crearon secuencias en la base de datos.")

        # Verificar que se generaron clusters en la base de datos
        clusters_count = self.clustering.session.query(Cluster).count()
        self.assertGreater(clusters_count, 0, "No se generaron clusters en la base de datos.")

        # Verificar que las entradas de los clusters están presentes
        cluster_entries_count = self.clustering.session.query(ClusterEntry).count()
        self.assertGreater(cluster_entries_count, 0, "No se generaron entradas de clusters.")

        # Verificar que las secuencias se agrupan correctamente en clusters
        clusters = self.clustering.session.query(Cluster).all()
        for cluster in clusters:
            entries = self.clustering.session.query(ClusterEntry).filter_by(cluster_id=cluster.id).count()
            self.assertGreater(entries, 0, f"El cluster {cluster.id} no tiene entradas asociadas.")

        # Verificar la coherencia de los datos almacenados
        cluster_entries = self.clustering.session.query(ClusterEntry).all()
        for entry in cluster_entries:
            sequence = self.clustering.session.get(Sequence, entry.sequence_id)
            self.assertIsNotNone(sequence, f"La secuencia con ID {entry.sequence_id} no existe.")

    if __name__ == '__main__':
        unittest.main()