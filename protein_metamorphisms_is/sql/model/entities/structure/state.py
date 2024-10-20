# state.py
from sqlalchemy import Column, Integer, String, ForeignKeyConstraint
from sqlalchemy.orm import relationship
from protein_metamorphisms_is.sql.model.core.base import Base

class State(Base):
    __tablename__ = 'state'

    id = Column(Integer, primary_key=True, nullable=False)
    file_path = Column(String, nullable=False)
    chain_id = Column(String, nullable=False)
    structure_id = Column(String, nullable=False)

    # Definir la clave foránea compuesta hacia Chain
    __table_args__ = (
        ForeignKeyConstraint(['chain_id', 'structure_id'], ['chain.id', 'chain.structure_id']),
    )

    # Relaciones
    chain = relationship("Chain", back_populates="states")

    def __repr__(self):
        return f"<State(id={self.id}, file_path={self.file_path}, chain_id={self.chain_id})>"
