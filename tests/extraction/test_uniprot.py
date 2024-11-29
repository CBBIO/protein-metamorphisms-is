import os
import unittest

import pytest
from sqlalchemy.dialects.mssql.information_schema import sequences

from protein_metamorphisms_is.helpers.config.yaml import read_yaml_config
from protein_metamorphisms_is.operation.extraction.uniprot import UniProtExtractor
from protein_metamorphisms_is.sql.model.entities.go_annotation.go_annotation import ProteinGOTermAnnotation
from protein_metamorphisms_is.sql.model.entities.go_annotation.go_term import GOTerm
from protein_metamorphisms_is.sql.model.entities.protein.accesion import Accession
from protein_metamorphisms_is.sql.model.entities.protein.protein import Protein
from protein_metamorphisms_is.sql.model.entities.sequence.sequence import Sequence
from protein_metamorphisms_is.sql.model.entities.structure.structure import Structure

@pytest.mark.order(2)
class TestUniProtExtractor(unittest.TestCase):

    def setUp(self):
        config_path = self.get_config_path()
        config = read_yaml_config(config_path)
        self.extractor = UniProtExtractor(config)

    def get_config_path(self):
        """Devuelve la ruta al archivo de configuración."""
        return os.path.join('tests/config/', "config.yaml")  # Ruta del archivo YAML de prueba

    def test_entities_created(self):
        """Verifica que se crean instancias en todas las entidades tras ejecutar el extractor."""
        # Ejecutar el extractor
        self.extractor.start()

        # Verificar GO Terms
        go_terms = self.extractor.session.query(GOTerm).all()
        self.assertGreater(len(go_terms), 0, "No se crearon términos GO en la base de datos.")

        # Verificar GO Annotations
        go_annotations = self.extractor.session.query(ProteinGOTermAnnotation).all()
        self.assertGreater(len(go_annotations), 0, "No se crearon anotaciones GO en la base de datos.")

        # Verificar Structures
        structures = self.extractor.session.query(Structure).all()
        self.assertGreater(len(structures), 0, "No se crearon estructuras en la base de datos.")

        # Verificar Sequences
        sequences = self.extractor.session.query(Sequence).all()
        self.assertGreater(len(sequences), 0, "No se crearon secuencias en la base de datos.")

        # Verificar Proteínas
        proteins = self.extractor.session.query(Protein).all()
        self.assertGreater(len(proteins), 0, "No se crearon proteínas en la base de datos.")

    if __name__ == '__main__':
        unittest.main()
