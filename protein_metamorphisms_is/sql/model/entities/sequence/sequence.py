from sqlalchemy import Column, Integer, String, Index, func
from sqlalchemy.orm import relationship
from protein_metamorphisms_is.sql.model.core.base import Base

# sequence.py
from sqlalchemy import Index, text

class Sequence(Base):
    __tablename__ = 'sequence'

    id = Column(Integer, primary_key=True, autoincrement=True)
    sequence = Column(String, nullable=False)
    sequence_hash = Column(String, unique=True, index=True)

    protein = relationship("Protein", back_populates="sequence", uselist=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if self.sequence and not self.sequence_hash:
            self.sequence_hash = func.md5(text(f"'{self.sequence}'"))

# Índice basado en sequence_hash para evitar el problema con sequence
Index('idx_sequence_hash', Sequence.sequence_hash)
