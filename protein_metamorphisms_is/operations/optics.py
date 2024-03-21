from sklearn.cluster import OPTICS
import numpy as np

from protein_metamorphisms_is.operations.base.operator import OperatorBase
from protein_metamorphisms_is.sql.model import Subcluster, SubclusterEntry, ClusterEntry, PDBChains, Cluster, \
    ChainEmbedding


class OpticsClustering(OperatorBase):
    def __init__(self, conf):
        super().__init__(conf)
        self.logger.info("OpticsClustering instance created")

    def start(self):
        try:
            self.logger.info("Starting OPTICS clustering process")
            cluster_ids = self.get_cluster_ids()
            for cluster_id in cluster_ids:
                embeddings, pdb_chain_ids = self.load_embeddings(cluster_id)
                if embeddings.size == 0:  # Verificar si no hay embeddings
                    continue
                cluster_labels = self.cluster_embeddings(embeddings)
                self.store_subclusters(cluster_id, cluster_labels, pdb_chain_ids)
            self.logger.info("Clustering process completed successfully")
        except Exception as e:
            self.logger.error(f"Error during clustering process: {e}")
            raise

    def get_cluster_ids(self):
        return [cluster.id for cluster in self.session.query(Cluster.id).all()]

    def load_embeddings(self, cluster_id):
        # Obtiene los embeddings y pdb_chain_ids para un cluster_id específico
        entries = self.session.query(
            ClusterEntry.pdb_chain_id,
            ChainEmbedding.embedding
        ).join(
            PDBChains, ClusterEntry.pdb_chain_id == PDBChains.id
        ).join(
            ChainEmbedding, PDBChains.id == ChainEmbedding.pdb_chain_id
        ).filter(
            ClusterEntry.cluster_id == cluster_id
        ).all()

        if not entries:
            return np.array([]), []

        # Asumiendo que 'embedding' es una lista o un numpy array ya compatible
        embeddings = np.array([entry.embedding for entry in entries])
        pdb_chain_ids = [entry.pdb_chain_id for entry in entries]

        return embeddings, pdb_chain_ids

    def cluster_embeddings(self, embeddings):
        # Ajusta min_samples basado en el número de muestras disponibles
        min_samples = min(5, len(embeddings) - 1)  # Asegura que min_samples nunca sea mayor que el número de muestras
        if min_samples < 2:  # OPTICS requiere al menos dos muestras para funcionar
            return np.array(
                [-1] * len(embeddings))  # Considera todos los puntos como ruido si no hay suficientes para clustering

        optics = OPTICS(min_samples=min_samples, xi=0.05, min_cluster_size=0.05)
        optics.fit(embeddings)
        return optics.labels_

    def store_subclusters(self, cluster_id, cluster_labels, pdb_chain_ids):
        # Diccionario para almacenar los subclusters y sus entradas
        subclusters_dict = {}

        # Iterar sobre cada etiqueta y pdb_chain_id juntos
        for label, pdb_chain_id in zip(cluster_labels, pdb_chain_ids):
            if label == -1:  # OPTICS puede marcar algunos puntos como ruido, los ignoramos
                continue

            if label not in subclusters_dict:
                subclusters_dict[label] = {
                    "entries": [],
                    "max_length": 0,
                    "representative_id": None
                }

            # Consulta la secuencia para el pdb_chain_id actual y calcula su longitud
            sequence = self.session.query(PDBChains.sequence).filter_by(id=pdb_chain_id).scalar()
            sequence_length = len(sequence)

            # Añadir la entrada al subcluster
            subclusters_dict[label]["entries"].append((pdb_chain_id, sequence_length))

            # Verificar si esta entrada tiene la secuencia de mayor longitud
            if sequence_length > subclusters_dict[label]["max_length"]:
                subclusters_dict[label]["max_length"] = sequence_length
                subclusters_dict[label]["representative_id"] = pdb_chain_id

        # Ahora, almacenar los subclusters y sus entradas en la base de datos
        for label, subcluster_info in subclusters_dict.items():
            subcluster = Subcluster(
                cluster_id=cluster_id,
                # Añadir otros campos necesarios aquí
            )
            self.session.add(subcluster)
            self.session.flush()  # Obtener el id del subcluster insertado

            # Marcar la entrada con la secuencia de mayor longitud como representativa
            for pdb_chain_id, sequence_length in subcluster_info["entries"]:
                is_representative = (pdb_chain_id == subcluster_info["representative_id"])
                subcluster_entry = SubclusterEntry(
                    subcluster_id=subcluster.id,
                    pdb_chain_id=pdb_chain_id,
                    is_representative=is_representative,
                    sequence_length=sequence_length,
                    # Añadir otros campos necesarios aquí
                )
                self.session.add(subcluster_entry)

        self.session.commit()



