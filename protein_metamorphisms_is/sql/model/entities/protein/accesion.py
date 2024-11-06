from sqlalchemy import Column, Integer, String, Boolean, DateTime, ForeignKey, func
from sqlalchemy.orm import relationship
from protein_metamorphisms_is.sql.model.core.base import Base

class Accession(Base):
    __tablename__ = "accession"

    id = Column(Integer, primary_key=True, autoincrement=True)
    accession_code = Column(String, unique=True, nullable=False)
    primary = Column(Boolean, default=True)
    tag = Column(String, nullable=True)
    protein_id = Column(String, ForeignKey('protein.id'))  # Relación con Protein
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())

    # Relaciones
    chain = relationship("Chain", back_populates="accession")
    protein = relationship("Protein", back_populates="accessions")
