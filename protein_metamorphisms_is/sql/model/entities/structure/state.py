from sqlalchemy import Column, Integer, String, ForeignKey
from sqlalchemy.orm import relationship
from protein_metamorphisms_is.sql.model.core.base import Base


class State(Base):
    __tablename__ = 'state'

    id = Column(Integer, primary_key=True, autoincrement=True, nullable=False)
    model_id = Column(String, nullable=False)
    file_path = Column(String, nullable=False)
    chain_id = Column(Integer, ForeignKey('chain.id'), nullable=False)  # Referencia al nuevo id de Chain
    structure_id = Column(String, nullable=False)

    # Relaciones
    chain = relationship("Chain", back_populates="states")

    def __repr__(self):
        return f"<State(id={self.id}, model_id={self.model_id}, file_path={self.file_path}, chain_id={self.chain_id})>"
