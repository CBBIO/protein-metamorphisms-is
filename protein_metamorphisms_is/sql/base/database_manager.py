from sqlalchemy import create_engine, QueuePool
from sqlalchemy.orm import sessionmaker

from protein_metamorphisms_is.sql.model.core.base import Base


class DatabaseManager:
    def __init__(self, conf):
        self.conf = conf
        self.engine = self.create_engine()
        self.Session = sessionmaker(bind=self.engine)

    def create_engine(self):
        """
        Create the SQLAlchemy engine with a QueuePool for connection management.
        """
        DATABASE_URI = (f"postgresql+psycopg2://{self.conf['DB_USERNAME']}:"
                        f"{self.conf['DB_PASSWORD']}"
                        f"@{self.conf['DB_HOST']}:"
                        f"{self.conf['DB_PORT']}/"
                        f"{self.conf['DB_NAME']}")
        # Crear el motor con el pool configurado
        engine = create_engine(
            DATABASE_URI,
            pool_size=0,          # Número máximo de conexiones en el pool
            max_overflow=0,       # Conexiones adicionales permitidas en caso de saturación
            poolclass=QueuePool,   # Clase de pool a utilizar
            pool_pre_ping=True     # Verificar conexiones antes de reutilizarlas
        )
        Base.metadata.create_all(engine)

        return engine

    def get_session(self):
        """
        Return a new session instance.
        """
        return self.Session()

    def get_engine(self):
        """
        Return the SQLAlchemy engine instance.
        """
        return self.engine

    def get_pool(self):
        """
        Return the pool object used by the engine.
        """
        return self.engine.pool
