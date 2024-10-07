import numpy as np
from sqlalchemy import text
from protein_metamorphisms_is.sql.model import SequenceEmbedding, PredictionMethod, \
    SequenceEmbeddingGOAnnotationTransfer, ClusterEntry
from protein_metamorphisms_is.tasks.queue import QueueTaskInitializer


class SequenceEmbeddingsGOAnnotationTransfer(QueueTaskInitializer):
    def __init__(self, conf):
        """
        Constructor de la clase que inicializa la configuración y el logger.
        """
        super().__init__(conf)
        self.logger.info("SequenceEmbeddingsGOAnnotationTransfer initialized")
        self.k = conf.get('k', 5)

    def enqueue(self):
        """
        Enqueue one task per cluster using the sequence of the representative.
        """
        try:
            # Obtener todos los clusters con su representante
            clusters = self.session.query(ClusterEntry).filter(ClusterEntry.is_representative == True).all()

            # Encolar una tarea por cada cluster con la secuencia del representante
            for cluster_entry in clusters:
                sequence = cluster_entry.sequence
                sequence_embedding = self.session.query(SequenceEmbedding).filter_by(sequence_id=sequence.id).first()

                # Si hay un embedding correspondiente para la secuencia del representante
                if sequence_embedding:
                    task_data = {
                        'sequence': sequence.sequence,
                        'sequence_id': sequence.id,
                        'embedding': sequence_embedding.embedding,
                        'cluster_id': cluster_entry.cluster_id
                    }
                    self.logger.info(
                        f"Enqueued task for cluster {cluster_entry.cluster_id} with representative sequence {sequence.id}")
                    self.publish_task(task_data)
                else:
                    self.logger.warning(f"No embedding found for sequence ID {sequence.id} in cluster {cluster_entry.cluster_id}")

        except Exception as e:
            self.logger.error(f"Error enqueuing tasks: {e}")

    def process(self, task_data):
        """
        Given the sequence and its corresponding embedding, find the K nearest clusters based on representative embeddings.
        """
        try:
            sequence = task_data['sequence']
            embedding = task_data['embedding']
            cluster_id = task_data['cluster_id']

            # Obtener todos los embeddings de los clusters diferentes al actual
            k_nearest_clusters = self.find_k_nearest_clusters(cluster_id, embedding)
            print(k_nearest_clusters)

            self.logger.info(f"Found {len(k_nearest_clusters)} nearest clusters for cluster ID {cluster_id}")
            return k_nearest_clusters

        except Exception as e:
            self.logger.error(f"Error processing task for cluster {task_data['cluster_id']}: {e}")
            return []

    import numpy as np
    from sqlalchemy import text

    import numpy as np
    from sqlalchemy import text

    def find_k_nearest_clusters(self, cluster_id, embedding):
        """
        Find the K nearest clusters based on the distance to the given embedding.
        Return the cluster IDs along with the calculated distance.
        """
        try:
            # Convertir el numpy array a lista antes de pasarlo como parámetro
            if isinstance(embedding, np.ndarray):
                embedding = embedding.tolist()

            query = text("""
                SELECT c.id AS nearest_cluster_id,
                       se.embedding <-> (
                           SELECT embedding
                           FROM sequence_embeddings
                           WHERE id = :cluster_entry_id
                       ) AS distance
                FROM clusters c
                JOIN cluster_entries ce ON c.id = ce.cluster_id
                JOIN sequence_embeddings se ON ce.sequence_id = se.sequence_id
                WHERE c.id != :cluster_id AND ce.is_representative = true
                ORDER BY distance
                LIMIT :k;
            """)

            result = self.session.execute(query, {
                'cluster_id': cluster_id,
                'cluster_entry_id': cluster_id,  # Pasa el ID del representante como parámetro
                'k': self.k
            }).fetchall()

            # Devuelve una lista de tuplas con (cluster_id, distancia)
            return [(r[0], r[1]) for r in result]

        except Exception as e:
            self.logger.error(f"Error finding nearest clusters for cluster ID {cluster_id}: {e}")
            return []

    def store_entry(self, predictions):
        """
        Store the prediction entries in the database.
        """
        try:
            print('Storing predictions...')
            print(predictions)
            # Implementación específica para almacenar los resultados
        except Exception as e:
            self.logger.error(f"Error storing predictions: {e}")
