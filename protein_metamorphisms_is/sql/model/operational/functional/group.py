from sqlalchemy import Column, Integer, DateTime, func, ForeignKey, String
from sqlalchemy.orm import relationship
from protein_metamorphisms_is.sql.model.core.base import Base


class GOTermPair(Base):
    __tablename__ = 'go_term_pair'

    id = Column(Integer, primary_key=True, autoincrement=True)
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())

    # Relaciones
    protein_id = Column(String, ForeignKey("protein.id"))  # Esta clave foránea vincula GOTermPair a Protein
    protein = relationship("Protein", back_populates="go_term_pairs")
    result = relationship("GOTermPairResult", back_populates="go_term_pair", uselist=False)
    entries = relationship("GOTermPairEntry", back_populates="go_term_pair", cascade="all, delete-orphan")

    def __repr__(self):
        return f"<GOTermPair(id={self.id})>"


class GOTermPairEntry(Base):
    __tablename__ = "go_term_pair_entry"

    id = Column(Integer, primary_key=True, autoincrement=True)
    go_term_pair_id = Column(Integer, ForeignKey("go_term_pair.id"), nullable=False)
    go_term_id = Column(String, ForeignKey("go_terms.go_id"), nullable=False)
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())

    # Relaciones
    go_term_pair = relationship("GOTermPair", back_populates="entries")
    go_term = relationship("GOTerm")

    def __repr__(self):
        return f"<GOTermPairEntry(go_term_pair_id={self.go_term_pair_id}, go_term_id={self.go_term_id})>"


class GOTermPairProtein(Base):
    __tablename__ = "go_term_pair_protein"

    id = Column(Integer, primary_key=True, autoincrement=True)
    go_term_pair_id = Column(Integer, ForeignKey("go_term_pair.id"), nullable=False)
    protein_id = Column(String, ForeignKey("protein.id"), nullable=False)
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())

    # Relaciones
    protein = relationship("Protein")
    go_term_pair = relationship("GOTermPair")

    def __repr__(self):
        return f"<GOTermPairProtein(go_term_pair_id={self.go_term_pair_id}, protein_id={self.protein_id})>"
