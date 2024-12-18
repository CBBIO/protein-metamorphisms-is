from sqlalchemy import Column, Integer, Float, ForeignKey, DateTime, func
from sqlalchemy.orm import relationship
from protein_metamorphisms_is.sql.model.core.base import Base


class GOTermPairResult(Base):
    __tablename__ = 'go_term_pair_result'

    id = Column(Integer, primary_key=True, autoincrement=True)
    go_term_pair_id = Column(Integer, ForeignKey('go_term_pair.id'), nullable=False, unique=True)
    mbl = Column(Float, nullable=True)
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())

    # Relaciones
    go_term_pair = relationship("GOTermPair", back_populates="result")

    def __repr__(self):
        return f"<GOTermPairResult(go_term_pair_id={self.go_term_pair_id})>"
