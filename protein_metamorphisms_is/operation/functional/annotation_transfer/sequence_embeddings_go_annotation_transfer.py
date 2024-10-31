from sqlalchemy import text

from protein_metamorphisms_is.sql.model.entities.embedding.sequence_embedding import SequenceEmbedding
from protein_metamorphisms_is.sql.model.entities.go_annotation.transference.cluster_go_term_annotation_transfer import \
    ClusterGOTermAnnotationTransfer
from protein_metamorphisms_is.sql.model.operational.clustering.cluster import ClusterEntry
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
            # Get all clusters with their representative
            clusters = self.session.query(ClusterEntry).filter(ClusterEntry.is_representative == True).all()

            # Enqueue a task for each cluster with the representative's sequence
            for cluster_entry in clusters[:50]:
                sequence = cluster_entry.sequence
                sequence_embeddings = self.session.query(SequenceEmbedding).filter_by(sequence_id=sequence.id).all()

                for sequence_embedding in sequence_embeddings:
                    if sequence_embedding:
                        task_data = {
                            'sequence': sequence.sequence,
                            'sequence_id': sequence.id,
                            'embedding': sequence_embedding.embedding,
                            'cluster_id': cluster_entry.cluster_id,
                            'embedding_type_id': sequence_embedding.embedding_type_id
                        }
                        self.logger.info(
                            f"Enqueued task for cluster {cluster_entry.cluster_id} with representative sequence {sequence.id}")
                        self.publish_task(task_data)
                    else:
                        self.logger.warning(
                            f"No embedding found for sequence ID {sequence.id} in cluster {cluster_entry.cluster_id}")


        except Exception as e:
            self.logger.error(f"Error enqueuing tasks: {e}")

    def process(self, task_data):
        try:
            cluster_id = task_data['cluster_id']
            embedding_type_id = task_data['embedding_type_id']
            task_type_id = task_data.get('task_type_id')

            query = text("""
            WITH TargetEmbedding AS (
    SELECT emb.embedding
    FROM "cluster" AS cl
    JOIN cluster_entry AS entry ON entry.cluster_id = cl.id
    JOIN "sequence" AS seq ON seq.id = entry.sequence_id
    JOIN sequence_embeddings AS emb ON emb.sequence_id = seq.id
    WHERE cl.id = :cluster_id 
      AND entry.is_representative = TRUE 
      AND emb.embedding_type_id = :embedding_type_id
),
RepresentativeEmbeddings AS (
    SELECT 
        cl.id AS cluster_id,
        emb.embedding
    FROM "cluster" AS cl
    JOIN cluster_entry AS entry ON entry.cluster_id = cl.id
    JOIN "sequence" AS seq ON seq.id = entry.sequence_id
    JOIN sequence_embeddings AS emb ON emb.sequence_id = seq.id
    WHERE emb.embedding_type_id = :embedding_type_id 
      AND entry.is_representative = TRUE
),
TopKClusters AS (
    SELECT 
        rep.cluster_id,
        l2_distance(target.embedding, rep.embedding) AS distance
    FROM 
        TargetEmbedding AS target, 
        RepresentativeEmbeddings AS rep
    ORDER BY 
        distance ASC
    LIMIT :k
),
GoTermAssociations AS (
    SELECT 
        entry.cluster_id,
        entry.id AS cluster_entry_id,
        top_k.distance,
        annotation.go_id AS go_term_id,
        prot.id AS protein_entry,
        ROW_NUMBER() OVER (PARTITION BY annotation.go_id ORDER BY top_k.distance ASC) AS rank
    FROM 
        TopKClusters AS top_k
    LEFT JOIN cluster_entry AS entry ON entry.cluster_id = top_k.cluster_id
    LEFT JOIN "sequence" AS seq ON seq.id = entry.sequence_id                
    LEFT JOIN "chain" AS ch ON ch.sequence_id = entry.sequence_id    
    LEFT JOIN "structure" AS struct ON struct.id = ch.structure_id   
    JOIN protein AS prot ON prot.id = struct.protein_id
    LEFT JOIN protein_go_term_annotation AS annotation ON annotation.protein_id = prot.id
    LEFT JOIN go_terms AS go ON go.go_id = annotation.go_id
)
SELECT 
    cluster_id,
    cluster_entry_id,
    distance,
    go_term_id,
    protein_entry
FROM 
    GoTermAssociations
WHERE 
    rank = 1  -- Seleccionar solo la fila con menor distancia por cada go_term_id
ORDER BY 
    cluster_id, cluster_entry_id;

            """)

            # Ejecutar la consulta
            with self.engine.connect() as connection:
                parameters = {'cluster_id': cluster_id, 'k': self.k, 'embedding_type_id': embedding_type_id}
                results = connection.execute(query, parameters).fetchall()
                self.logger.info(f"Query executed successfully for cluster_id: {cluster_id}")

            # Filtrar y almacenar predicciones
            predictions = []
            for result in results:
                prediction = {
                    'go_term_id': result.go_term_id,
                    'source_cluster_id': cluster_id,
                    'target_cluster_id': result.cluster_id,
                    'distance': result.distance,
                    'embedding_type_id': embedding_type_id,
                    'protein_entry': result.protein_entry,
                    'task_type_id': task_type_id
                }
                predictions.append(prediction)

            return predictions

        except Exception as e:
            self.logger.error(f"Unexpected error in process: {e}")
            return []

    def store_entry(self, predictions):
        """
        Almacena las entradas de predicción en la base de datos, verificando si ya existen.
        """
        try:
            self.logger.info("Storing predictions...")

            for prediction in predictions:
                # Validar que 'go_id' no sea None para evitar violación de restricciones
                if prediction.get('go_term_id') is None:
                    self.logger.warning(f"Skipping entry due to missing go_id: {prediction}")
                    continue

                # Crear una nueva instancia de ClusterGOTermAnnotationTransfer con los datos de predicción
                transfer_entry = ClusterGOTermAnnotationTransfer(
                    go_id=prediction['go_term_id'],
                    protein_entry_name=prediction['protein_entry'],
                    source_cluster_id=prediction['source_cluster_id'],
                    target_cluster_id=prediction['target_cluster_id'],
                    distance=prediction['distance'],
                    is_transferred=prediction.get('is_transferred', True),
                    embedding_type_id=prediction['embedding_type_id']
                )

                # Añadir la entrada a la sesión
                self.session.add(transfer_entry)
                self.logger.info(
                    f"Stored new entry for protein: {prediction['protein_entry']}, GO term: {prediction['go_term_id']}")

            # Confirmar la transacción para guardar los registros en la base de datos
            self.session.commit()
            self.logger.info("Predictions successfully stored.")

        except Exception as e:
            self.logger.error(f"Error storing predictions: {e}")
            self.session.rollback()



