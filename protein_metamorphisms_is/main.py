import argparse
from protein_metamorphisms_is.helpers.config.yaml import read_yaml_config

from protein_metamorphisms_is.sql.model.model import (
    AccessionManager,
    UniProtExtractor,
    SequenceEmbeddingManager
)


def main(config_path='config/config.yaml'):
    conf = read_yaml_config(config_path)

    AccessionManager(conf).fetch_accessions_from_api()
    # AccessionManager(conf).load_accessions_from_csv()
    UniProtExtractor(conf).start()
    SequenceEmbeddingManager(conf).start()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the pipeline with a specified configuration file.")
    parser.add_argument("--config", type=str, required=False, help="Path to the configuration YAML file.")
    args = parser.parse_args()

    config_path = args.config if args.config else 'config/config.yaml'
    main(config_path=config_path)
