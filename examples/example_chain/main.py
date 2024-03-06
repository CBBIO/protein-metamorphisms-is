from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from protein_data_handler.chain import ChainExtractor
from protein_data_handler.helpers.config.yaml import read_yaml_config
from protein_data_handler.sql.model import Base

if __name__ == "__main__":
    config = read_yaml_config("./config.yaml")
    DATABASE_URI = \
        (f"postgresql+psycopg2://{config['DB_USERNAME']}:"
         f"{config['DB_PASSWORD']}"
         f"@{config['DB_HOST']}:"
         f"{config['DB_PORT']}/"
         f"{config['DB_NAME']}")
    engine = create_engine(DATABASE_URI)
    Session = sessionmaker(bind=engine)
    session = Session()
    Base.metadata.create_all(engine)
    chain_extractor = ChainExtractor(session,config['pdb_dir'],config['chain_dir'])
    chain_extractor.create_decomposed_structure()

