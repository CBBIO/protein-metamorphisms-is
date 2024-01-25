from protein_data_handler.helpers.config.yaml import read_yaml_config
from protein_data_handler.information_system.pdb import PDBExtractor
from protein_data_handler.information_system.uniprot import UniProtExtractor
from protein_data_handler.operations.cdhit import CDHit
from protein_data_handler.operations.cealign import CEAlign

if __name__ == "__main__":
    conf = read_yaml_config("./config.yaml")

    # UniProtExtractor(conf).start()

    # PDBExtractor(conf).start()

    # CDHit(conf).start()

    CEAlign(conf).start()
