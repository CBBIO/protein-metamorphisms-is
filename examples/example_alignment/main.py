from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from protein_data_handler.alignment import UniProtPDBMapping
from protein_data_handler.helpers.config.yaml import read_yaml_config

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

    mapping = UniProtPDBMapping(session)
    pares = mapping.realizar_consulta_cadenas_iguales()
    mapping.volcar_datos_alineamiento(pares)

