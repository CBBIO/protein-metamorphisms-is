import os
import unittest

import pytest

from protein_metamorphisms_is.helpers.config.yaml import read_yaml_config
from protein_metamorphisms_is.operation.extraction.pdb import PDBExtractor
from protein_metamorphisms_is.sql.model.entities.protein.accesion import Accession
from protein_metamorphisms_is.sql.model.entities.protein.protein import Protein
from protein_metamorphisms_is.sql.model.entities.sequence.sequence import Sequence
from protein_metamorphisms_is.sql.model.entities.structure.structure import Structure
from protein_metamorphisms_is.sql.model.entities.structure.chain import Chain
from protein_metamorphisms_is.sql.model.entities.structure.state import State

@pytest.mark.order(3)
class TestPDBExtractor(unittest.TestCase):

    def setUp(self):
        config_path = self.get_config_path()
        config = read_yaml_config(config_path)
        self.extractor = PDBExtractor(config)

    def get_config_path(self):
        """Devuelve la ruta al archivo de configuración."""
        return os.path.join('tests/config/', "config.yaml")  # Ruta del archivo YAML de prueba

    def test_entities_created(self):
        """Verifica que se crean instancias en todas las entidades tras ejecutar el extractor."""
        # Ejecutar el extractor
        self.extractor.start()

        # Verificar Structures
        structures = self.extractor.session.query(Structure).all()
        self.assertGreater(len(structures), 0, "No se crearon estructuras en la base de datos.")

        # Verificar Chains
        chains = self.extractor.session.query(Chain).all()
        self.assertGreater(len(chains), 0, "No se crearon cadenas en la base de datos.")

        # Verificar States
        states = self.extractor.session.query(State).all()
        self.assertGreater(len(states), 0, "No se crearon estados en la base de datos.")

        # Verificar Sequences
        sequences = self.extractor.session.query(Sequence).all()
        self.assertGreater(len(sequences), 0, "No se crearon secuencias en la base de datos.")

        # Verificar Proteínas
        proteins = self.extractor.session.query(Protein).all()
        self.assertGreater(len(proteins), 0, "No se crearon proteínas en la base de datos.")

        # Verificar Accessions
        accessions = self.extractor.session.query(Accession).all()
        self.assertGreater(len(accessions), 0, "No se crearon accesiones en la base de datos.")

    if __name__ == '__main__':
        unittest.main()
