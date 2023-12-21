from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from protein_data_handler.helpers.config.yaml import read_yaml_config
from protein_data_handler.models.uniprot import Base
from protein_data_handler.structural_alignment import CEAlignHandler

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
    ce_align = CEAlignHandler(session)
    ce_align.get_targets()
    ce_align.all_to_all_align(identity_threshold=config['identity_threshold'])
