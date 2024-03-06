from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from protein_metamorphisms_is.sql.model import Base


class DatabaseManager:
    def __init__(self, conf):
        self.conf = conf
        self.engine = self.create_engine()
        self.Session = sessionmaker(bind=self.engine)

    def create_engine(self):
        DATABASE_URI = (f"postgresql+psycopg2://{self.conf['DB_USERNAME']}:"
                        f"{self.conf['DB_PASSWORD']}"
                        f"@{self.conf['DB_HOST']}:"
                        f"{self.conf['DB_PORT']}/"
                        f"{self.conf['DB_NAME']}")
        engine = create_engine(DATABASE_URI, pool_size=0, max_overflow=0)
        Base.metadata.create_all(engine)

        return engine

    def get_session(self):
        return self.Session()

    def get_engine(self):
        return self.engine
