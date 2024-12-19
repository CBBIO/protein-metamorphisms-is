from datetime import datetime

from sqlalchemy import Column, Integer, DateTime, ForeignKey, Boolean, Float, String
from sqlalchemy.orm import relationship

from protein_metamorphisms_is.sql.model.core.base import Base


class Cluster(Base):
    __tablename__ = 'cluster'
    id = Column(Integer, primary_key=True)
    created_at = Column(DateTime, default=datetime.now)

    # Relación con ClusterEntries
    entries = relationship("ClusterEntry", back_populates="cluster", cascade="all, delete-orphan")

    # Relación con Subclusters
    subclusters = relationship("Subcluster", back_populates="cluster", cascade="all, delete-orphan")

    def __repr__(self):
        return f"<Cluster(id={self.id})>"


class ClusterEntry(Base):
    __tablename__ = 'cluster_entry'
    id = Column(Integer, primary_key=True)
    cluster_id = Column(Integer, ForeignKey('cluster.id', ondelete="CASCADE"), nullable=False)
    sequence_id = Column(Integer, ForeignKey('sequence.id'), nullable=False)
    is_representative = Column(Boolean)
    sequence_length = Column(Integer)
    identity = Column(Float)
    created_at = Column(DateTime, default=datetime.now)

    # Relaciones
    cluster = relationship("Cluster", back_populates="entries")
    sequence = relationship("Sequence")


class Subcluster(Base):
    __tablename__ = 'subcluster'
    id = Column(Integer, primary_key=True)
    cluster_id = Column(Integer, ForeignKey('cluster.id', ondelete="CASCADE"), nullable=False)
    description = Column(String, nullable=True)
    created_at = Column(DateTime, default=datetime.now)

    # Relaciones
    cluster = relationship("Cluster", back_populates="subclusters")
    subcluster_entries = relationship("SubclusterEntry", back_populates="subcluster", cascade="all, delete-orphan")


class SubclusterEntry(Base):
    __tablename__ = 'subcluster_entry'
    id = Column(Integer, primary_key=True)
    subcluster_id = Column(Integer, ForeignKey('subcluster.id', ondelete="CASCADE"), nullable=False)
    structure_3di_id = Column(Integer, ForeignKey('structure_3di.id'), nullable=False)
    is_representative = Column(Boolean, default=False)
    created_at = Column(DateTime, default=datetime.now)
    sequence_length = Column(Integer)
    identity = Column(Float)

    # Relaciones
    subcluster = relationship("Subcluster", back_populates="subcluster_entries")
    structure_3di = relationship("Structure3Di")
    alignment_group_entry = relationship("AlignmentGroupEntry", uselist=False, back_populates="subcluster_entry")  # Relación uno a uno
