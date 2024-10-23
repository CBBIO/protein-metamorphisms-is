from sqlalchemy import Column, String, Integer, ForeignKey
from sqlalchemy.orm import relationship
from protein_metamorphisms_is.sql.model.core.base import Base

class Chain(Base):
    __tablename__ = 'chain'

    id = Column(Integer, primary_key=True, autoincrement=True)  # Autoincremental
    name = Column(String, nullable=False)  # Nombre de la cadena
    structure_id = Column(String, ForeignKey('structure.id'), nullable=False)
    sequence_id = Column(Integer, ForeignKey('sequence.id'))

    # Relaciones
    structure = relationship("Structure", back_populates="chains")
    states = relationship("State", back_populates="chain", cascade="all, delete-orphan")
    sequence = relationship("Sequence", back_populates="chain")

    def __repr__(self):
        return f"<Chain(id={self.id}, name={self.name}, structure_id={self.structure_id})>"
