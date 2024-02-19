from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker


def create_session(database_uri):
    """
    Crea y devuelve una sesión de SQLAlchemy.

    :param database_uri: URI de la base de datos para conectar.
    :return: Instancia de sesión de SQLAlchemy.
    """

    engine = create_engine(database_uri)
    Session = sessionmaker(bind=engine)
    return Session()
