from sqlalchemy import Column, Integer, String, Float, Boolean, ForeignKey
from sqlalchemy.orm import relationship

from protein_metamorphisms_is.sql.model.core.base import Base


class ClusterGOTermAnnotationTransfer(Base):
    __tablename__ = 'cluster_go_term_annotation_transfer'

    id = Column(Integer, primary_key=True)
    go_id = Column(String, nullable=False)
    protein_entry_name = Column(String, nullable=False)
    source_cluster_id = Column(Integer, ForeignKey('cluster.id'))
    target_cluster_id = Column(Integer, ForeignKey('cluster.id'))
    distance = Column(Float, nullable=False)
    is_transferred = Column(Boolean, default=True)
    embedding_type_id = Column(Integer, nullable=False)

    # Optional relationships if needed
    source_cluster = relationship("Cluster", foreign_keys=[source_cluster_id])
    target_cluster = relationship("Cluster", foreign_keys=[target_cluster_id])
