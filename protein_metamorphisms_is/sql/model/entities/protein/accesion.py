from sqlalchemy import Column, Integer, String, Boolean, DateTime, func

from protein_metamorphisms_is.sql.model.core.base import Base


class Accession(Base):
    """
    Represents a basic accession entity for tracking proteins.
    """

    __tablename__ = "accession"

    id = Column(Integer, primary_key=True, autoincrement=True)
    accession_code = Column(String, unique=True, nullable=False)
    primary = Column(Boolean, default=True)
    tag = Column(String, nullable=True)
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())
