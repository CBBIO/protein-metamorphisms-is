from sklearn.cluster import OPTICS
import numpy as np

from protein_metamorphisms_is.operation.base.operator import OperatorBase
from protein_metamorphisms_is.sql.model import Subcluster, SubclusterEntry, ClusterEntry, SequenceEmbedding, Sequence, EmbeddingType, Cluster

class OpticsClustering(OperatorBase):
    def __init__(self, conf):
        super().__init__(conf)
        self.logger.info("OpticsClustering instance created")

    def start(self):
        try:
            self.logger.info("Starting OPTICS clustering process")
            cluster_ids = self.get_cluster_ids()
            embedding_types = self.session.query(EmbeddingType).all()
            for cluster_id in cluster_ids:
                for embedding_type in embedding_types:
                    embeddings, sequence_ids = self.load_embeddings(cluster_id, embedding_type)
                    if embeddings.size == 0:
                        continue
                    cluster_labels = self.cluster_embeddings(embeddings)
                    print(cluster_labels)
                    self.store_subclusters(cluster_id, cluster_labels, sequence_ids, embedding_type)
            self.logger.info("Clustering process completed successfully")
        except Exception as e:
            self.logger.error(f"Error during clustering process: {e}")
            raise

    def get_cluster_ids(self):
        return [cluster.id for cluster in self.session.query(Cluster.id).all()]

    def load_embeddings(self, cluster_id, embedding_type):
        # Get embeddings and sequence_ids for a specific cluster_id
        entries = self.session.query(
            ClusterEntry.sequence_id,
            SequenceEmbedding.embedding
        ).join(
            SequenceEmbedding, ClusterEntry.sequence_id == SequenceEmbedding.sequence_id
        ).filter(
            ClusterEntry.cluster_id == cluster_id,
            SequenceEmbedding.embedding_type == embedding_type
        ).all()

        if not entries:
            return np.array([]), []
        embeddings = np.array([entry.embedding for entry in entries])
        sequence_ids = [entry.sequence_id for entry in entries]
        print('entra')
        return embeddings, sequence_ids

    def cluster_embeddings(self, embeddings):
        # Adjust min_samples based on the number of samples available
        min_samples = min(5, len(embeddings) - 1)
        if min_samples < 2:
            return np.array([-1] * len(embeddings))

        optics = OPTICS(min_samples=min_samples, xi=0.000001, min_cluster_size=0.05)
        optics.fit(embeddings)
        print(optics.labels_)
        return optics.labels_

    def store_subclusters(self, cluster_id, cluster_labels, sequence_ids, embedding_type):
        subclusters_dict = {}

        for label, sequence_id in zip(cluster_labels, sequence_ids):
            if label == -1:
                continue

            if label not in subclusters_dict:
                subclusters_dict[label] = {
                    "entries": [],
                    "max_length": 0,
                    "representative_id": None
                }

            sequence = self.session.query(Sequence.sequence).filter(Sequence.id == sequence_id).first()
            sequence_length = len(sequence[0])

            subclusters_dict[label]["entries"].append((sequence_id, sequence_length))

            if sequence_length > subclusters_dict[label]["max_length"]:
                subclusters_dict[label]["max_length"] = sequence_length
                subclusters_dict[label]["representative_id"] = sequence_id

        for label, subcluster_info in subclusters_dict.items():
            subcluster = Subcluster(
                cluster_id=cluster_id,
                embedding_type=embedding_type
            )
            self.session.add(subcluster)
            self.session.flush()

            for sequence_id, sequence_length in subcluster_info["entries"]:
                is_representative = (sequence_id == subcluster_info["representative_id"])
                subcluster_entry = SubclusterEntry(
                    subcluster_id=subcluster.id,
                    sequence_id=sequence_id,
                    is_representative=is_representative,
                    sequence_length=sequence_length,
                )
                self.session.add(subcluster_entry)

        self.session.commit()

