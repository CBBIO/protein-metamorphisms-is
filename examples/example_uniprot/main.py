from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from protein_data_handler.helpers.config.yaml import read_yaml_config
from protein_data_handler.models.uniprot import Base
from protein_data_handler.uniprot import cargar_codigos_acceso, extraer_entradas
import logging
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


if __name__ == "__main__":
    config = read_yaml_config("./config.yaml")
    DATABASE_URI = \
        (f"postgresql+psycopg2://{config['DB_USERNAME']}:"
         f"{config['DB_PASSWORD']}"
         f"@{config['DB_HOST']}:"
         f"{config['DB_PORT']}/"
         f"{config['DB_NAME']}")
    engine = create_engine(DATABASE_URI)
    Base.metadata.drop_all(engine)

    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()

    cargar_codigos_acceso(
        criterio_busqueda=config['criterio_busqueda'], limite=config['limit'], session=session)
    extraer_entradas(session=session)
