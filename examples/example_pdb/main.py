from protein_data_handler.helpers.config.yaml import read_yaml_config
from protein_data_handler.pdb import download_entire_pdb
import logging

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

config = read_yaml_config("./config.yaml")

server = config['server']
pdb = config['pdb']
file_format = config['file_format']


if __name__ == "__main__":
    download_entire_pdb(server, pdb, file_format)
