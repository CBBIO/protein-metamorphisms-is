from sqlalchemy import Column, String, DateTime, func, ForeignKey, Float
from sqlalchemy.orm import relationship
from protein_metamorphisms_is.sql.model.core.base import Base


class Structure(Base):
    __tablename__ = 'structure'

    id = Column(String, primary_key=True, unique=True, nullable=False)
    protein_id = Column(String, ForeignKey('protein.id'), nullable=False)
    method = Column(String)
    resolution = Column(Float)
    file_path = Column(String, nullable=False, unique=True)
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())

    # Relationships
    protein = relationship("Protein", back_populates="structure")
    chains = relationship("Chain", back_populates="structure", cascade="all, delete-orphan")  # One-to-many relationship

    def __repr__(self):
        return f"<Structure(id={self.id}, protein_id={self.protein_id}, method={self.method}, resolution={self.resolution})>"
