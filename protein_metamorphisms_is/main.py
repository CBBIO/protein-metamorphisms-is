import logging


def main(config_path='config/config.yaml'):
    from protein_metamorphisms_is.helpers.config.yaml import read_yaml_config
    conf = read_yaml_config(config_path)

    # ✅ Mover esto al principio
    logger = logging.getLogger("protein_metamorphisms_is")
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler()
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    # Step 1: Ensure pgvector is available before ORM definitions
    from protein_metamorphisms_is.sql.base.database_manager import DatabaseManager
    DatabaseManager(conf)

    # Step 2: Import ORM-based logic after pgvector is ready
    from protein_metamorphisms_is.sql.model.model import (
        AccessionManager,
        UniProtExtractor,
        SequenceEmbeddingManager
    )

    # Step 3: Run components
    AccessionManager(conf).load_accessions_from_csv()
    UniProtExtractor(conf).start()
    SequenceEmbeddingManager(conf).start()


if __name__ == '__main__':
    main()
