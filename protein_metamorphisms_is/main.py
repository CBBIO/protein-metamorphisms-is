from protein_metamorphisms_is.helpers.config.yaml import read_yaml_config
from protein_metamorphisms_is.information_system.pdb import PDBExtractor
from protein_metamorphisms_is.information_system.uniprot import UniProtExtractor
from protein_metamorphisms_is.operations.cdhit import CDHit
from protein_metamorphisms_is.operations.structural_alignment import StructuralAlignmentManager


if __name__ == "__main__":
    conf = read_yaml_config("config/config.yaml")

    UniProtExtractor(conf).start()

    PDBExtractor(conf).start()

    CDHit(conf).start()

    StructuralAlignmentManager(conf).start()


