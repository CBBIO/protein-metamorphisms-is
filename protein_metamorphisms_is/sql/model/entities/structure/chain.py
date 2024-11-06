from sqlalchemy import Column, String, Integer, ForeignKey
from sqlalchemy.orm import relationship
from protein_metamorphisms_is.sql.model.core.base import Base

class Chain(Base):
    __tablename__ = 'chain'

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String, nullable=False)
    structure_id = Column(String, ForeignKey('structure.id'), nullable=False)
    sequence_id = Column(Integer, ForeignKey('sequence.id'))
    accession_id = Column(Integer, ForeignKey('accession.id'), unique=False)  # Relación uno a uno con Accession

    # Relaciones
    structure = relationship("Structure", back_populates="chains")
    states = relationship("State", back_populates="chain", cascade="all, delete-orphan")
    sequence = relationship("Sequence", back_populates="chain")
    accession = relationship("Accession", back_populates="chain", uselist=False)  # Relación uno a uno con Accession
