from abc import abstractmethod, ABC

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from protein_data_handler.helpers.logger.logger import setup_logger
from protein_data_handler.sql.model import Base


class BioinfoExtractorBase(ABC):
    """
    An abstract base class for operating bioinformatics data.

    This class provides a framework for connecting to and interacting with various
    bioinformatics data sources. It is designed to be subclassed with specific
    implementations for different data sources.

    Attributes:
        conf (dict): Configuration parameters for the extractor.
        logger (Logger): Logger object for logging information.
        session (Session, optional): A SQLAlchemy session for database interactions.

    Args:
        conf (dict): Configuration dictionary containing necessary parameters.
        session_required (bool): Flag to indicate if a database session is required.
    """

    def __init__(self, conf, session_required):
        """
        Initialize the BioinfoExtractorBase class.

        Sets up the configuration and logger. Initializes the database session if required.
        """
        self.conf = conf
        self.logger = setup_logger(self.__class__.__name__)
        self.logger.info(f"Initializing {self.__class__.__name__}")

        if session_required:
            self.session_init()

    def session_init(self):
        """
        Initialize the database session.

        Sets up the database connection using SQLAlchemy and creates the database schema.
        """
        self.logger.info("Initializing database session")
        DATABASE_URI = \
            (f"postgresql+psycopg2://{self.conf['DB_USERNAME']}:"
             f"{self.conf['DB_PASSWORD']}"
             f"@{self.conf['DB_HOST']}:"
             f"{self.conf['DB_PORT']}/"
             f"{self.conf['DB_NAME']}")
        engine = create_engine(DATABASE_URI)
        self.engine = engine
        Base.metadata.drop_all(engine)
        Base.metadata.create_all(engine)
        Session = sessionmaker(bind=engine)

        self.session = Session()

    @abstractmethod
    def start(self):
        """
        Start the data extraction process.

        This abstract method should be implemented by all subclasses to define
        the specific data extraction logic for each bioinformatics data source.
        """
        pass
