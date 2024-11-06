from sqlalchemy import text

from protein_metamorphisms_is.sql.model.entities.embedding.sequence_embedding import SequenceEmbedding
from protein_metamorphisms_is.sql.model.entities.go_annotation.transference.sequence_go_term_annotation import \
    SequenceGoTermAnnotation
from protein_metamorphisms_is.sql.model.operational.clustering.cluster import ClusterEntry
from protein_metamorphisms_is.tasks.queue import QueueTaskInitializer




class SequenceGOAnnotation(QueueTaskInitializer):
    def __init__(self, conf):
        """
        Constructor de la clase que inicializa la configuración y el logger.
        """
        super().__init__(conf)
        self.logger.info("SequenceEmbeddingsGOAnnotationTransfer initialized")
        self.k = conf.get('k', 5)

    def enqueue(self):
        """
        Enqueue a task for each sequence instead of clusters.
        """
        try:
            # Get all sequences
            sequences = self.session.query(SequenceEmbedding).all()

            for sequence_embedding in sequences:
                task_data = {
                    'sequence': sequence_embedding.sequence.sequence,
                    'sequence_id': sequence_embedding.sequence_id,
                    'embedding': sequence_embedding.embedding,
                    'embedding_type_id': sequence_embedding.embedding_type_id
                }
                self.logger.info(f"Enqueued task for sequence {sequence_embedding.sequence_id}")
                self.publish_task(task_data)

        except Exception as e:
            self.logger.error(f"Error enqueuing tasks: {e}")

    def process(self, task_data):
        try:
            sequence_id = task_data['sequence_id']
            embedding_type_id = task_data['embedding_type_id']
            task_type_id = task_data.get('task_type_id')

            query = text("""
            WITH TargetEmbedding AS (
                SELECT embedding
                FROM sequence_embeddings
                WHERE sequence_id = :sequence_id AND embedding_type_id = :embedding_type_id
            ),
            AllSequenceEmbeddings AS (
                SELECT 
                    seq.id AS sequence_id,
                    seq_embedding.embedding
                FROM "sequence" AS seq
                JOIN sequence_embeddings AS seq_embedding ON seq_embedding.sequence_id = seq.id
                WHERE seq_embedding.embedding_type_id = :embedding_type_id
            ),
            TopKSequences AS (
                SELECT 
                    all_seq.sequence_id,
                    l2_distance(target.embedding, all_seq.embedding) AS distance
                FROM 
                    TargetEmbedding AS target, 
                    AllSequenceEmbeddings AS all_seq
                ORDER BY 
                    distance ASC
                LIMIT :k
            ),
            GoTermAssociations AS (
                SELECT 
                    seq.id AS sequence_id,
                    top_k.distance,
                    annotation.go_id AS go_term_id,
                    prot.id AS protein_entry,
                    ROW_NUMBER() OVER (PARTITION BY annotation.go_id ORDER BY top_k.distance ASC) AS rank
                FROM 
                    TopKSequences AS top_k
                LEFT JOIN "sequence" AS seq ON seq.id = top_k.sequence_id
                LEFT JOIN "chain" AS ch ON ch.sequence_id = seq.id    
                LEFT JOIN "structure" AS struct ON struct.id = ch.structure_id   
                JOIN protein AS prot ON prot.id = struct.protein_id
                LEFT JOIN protein_go_term_annotation AS annotation ON annotation.protein_id = prot.id
                LEFT JOIN go_terms AS go ON go.go_id = annotation.go_id
            )
            SELECT 
                sequence_id,
                distance,
                go_term_id,
                protein_entry
            FROM 
                GoTermAssociations
            WHERE 
                rank = 1  -- Seleccionar solo la fila con menor distancia por cada go_term_id
            ORDER BY 
                sequence_id;
            """)

            # Ejecutar la consulta
            with self.engine.connect() as connection:
                parameters = {'sequence_id': sequence_id, 'k': self.k, 'embedding_type_id': embedding_type_id}
                results = connection.execute(query, parameters).fetchall()
                self.logger.info(f"Query executed successfully for sequence_id: {sequence_id}")

            # Filtrar y almacenar predicciones
            predictions = []
            for result in results:
                prediction = {
                    'go_term_id': result.go_term_id,
                    'source_sequence_id': sequence_id,
                    'target_sequence_id': result.sequence_id,
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
            return
            self.logger.info("Storing predictions...")

            for prediction in predictions:
                if prediction.get('go_term_id') is None:
                    self.logger.warning(f"Skipping entry due to missing go_id: {prediction}")
                    continue

                transfer_entry = ProteinGOTermAnnotation(
                    go_id=prediction['go_term_id'],
                    protein_entry_name=prediction['protein_entry'],
                    source_sequence_id=prediction['source_sequence_id'],
                    target_sequence_id=prediction['target_sequence_id'],
                    distance=prediction['distance'],
                    is_transferred=prediction.get('is_transferred', True),
                    embedding_type_id=prediction['embedding_type_id']
                )

                self.session.add(transfer_entry)
                self.logger.info(
                    f"Stored new entry for protein: {prediction['protein_entry']}, GO term: {prediction['go_term_id']}")

            self.session.commit()
            self.logger.info("Predictions successfully stored.")

        except Exception as e:
            self.logger.error(f"Error storing predictions: {e}")
            self.session.rollback()


