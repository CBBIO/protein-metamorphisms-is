# chain.py
from sqlalchemy import Column, String, ForeignKey
from sqlalchemy.orm import relationship
from protein_metamorphisms_is.sql.model.core.base import Base
from protein_metamorphisms_is.sql.model.entities.structure.state import State

class Chain(Base):
    __tablename__ = 'chain'

    id = Column(String, primary_key=True)  # Nombre de la cadena
    structure_id = Column(String, ForeignKey('structure.id'), primary_key=True)  # ID de la estructura

    # Relaciones
    structure = relationship("Structure", back_populates="chains")
    states = relationship("State", back_populates="chain", cascade="all, delete-orphan")

    def __repr__(self):
        return f"<Chain(id={self.id}, structure_id={self.structure_id})>"
