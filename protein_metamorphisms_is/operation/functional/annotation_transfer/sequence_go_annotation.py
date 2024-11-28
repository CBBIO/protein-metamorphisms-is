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
            WITH target_embedding AS (
                SELECT embedding
                FROM sequence_embeddings
                WHERE sequence_id = :sequence_id AND embedding_type_id = :embedding_type_id
            ),
            filtered_sequences AS (
                SELECT 
                    se.sequence_id,
                    se.embedding,
                    (se.embedding <-> te.embedding) AS distance
                FROM 
                    sequence_embeddings se,
                    target_embedding te
                WHERE 
                    (se.embedding <-> te.embedding) < 1
            ),
            ranked_sequences AS (
                SELECT 
                    fs.sequence_id,
                    fs.distance,
                    COALESCE(p.id, acc.protein_id) AS protein_id,
                    COALESCE(p.description, 'Derived from Chain') AS protein_description,
                    c.id AS chain_id,
                    c.name AS chain_name
                FROM 
                    filtered_sequences fs
                LEFT JOIN 
                    chain c ON c.sequence_id = fs.sequence_id
                LEFT JOIN 
                    accession acc ON acc.code = c.accession_code
                LEFT JOIN 
                    protein p ON p.sequence_id = fs.sequence_id OR p.id = acc.protein_id
                ORDER BY 
                    fs.distance ASC
            ),
            ranked_go_terms AS (
                SELECT 
                    rs.sequence_id,
                    rs.distance,
                    rs.protein_id,
                    rs.protein_description,
                    rs.chain_id,
                    rs.chain_name,
                    gta.go_id,
                    gt.category AS go_category,
                    gt.description AS go_description,
                    ROW_NUMBER() OVER (PARTITION BY gta.go_id ORDER BY rs.distance ASC) AS rank
                FROM ranked_sequences rs
                LEFT JOIN protein_go_term_annotation gta ON gta.protein_id = rs.protein_id
                LEFT JOIN go_terms gt ON gt.go_id = gta.go_id
            )
            SELECT 
                sequence_id,
                distance,
                protein_id,
                protein_description,
                chain_id,
                chain_name,
                go_id,
                go_category,
                go_description
            FROM ranked_go_terms
            WHERE rank = 1
            ORDER BY distance ASC, protein_id;
            """)

            # Ejecutar la consulta
            with self.engine.connect() as connection:
                parameters = {'sequence_id': sequence_id, 'embedding_type_id': embedding_type_id}
                results = connection.execute(query, parameters).fetchall()
                self.logger.info(f"Query executed successfully for sequence_id: {sequence_id}")

            # Filtrar y almacenar predicciones
            predictions = []
            for result in results:
                prediction = {
                    'go_term_id': result.go_id,
                    'source_sequence_id': sequence_id,
                    'target_sequence_id': result.sequence_id,
                    'distance': result.distance,
                    'embedding_type_id': embedding_type_id,
                    'protein_entry': result.protein_id,
                    'task_type_id': task_type_id
                }
                predictions.append(prediction)
            print(predictions)
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
                if not prediction.get('go_term_id'):
                    self.logger.warning(f"Skipping entry due to missing go_id: {prediction}")
                    continue

                # Verificar si la entrada ya existe
                existing_entry = self.session.query(SequenceGoTermAnnotation).filter_by(
                    go_id=prediction['go_term_id'],
                    sequence_id=prediction['source_sequence_id']  # Relación con la secuencia correspondiente
                ).first()

                if existing_entry:
                    self.logger.info(
                        f"Entry already exists for GO term: {prediction['go_term_id']} and sequence: {prediction['target_sequence_id']}")
                    continue

                # Crear nueva entrada
                transfer_entry = SequenceGoTermAnnotation(
                    go_id=prediction['go_term_id'],
                    sequence_id=prediction['source_sequence_id'],  # Relación con la secuencia
                    distance=prediction['distance']
                )

                self.session.add(transfer_entry)
                self.logger.info(
                    f"Stored new entry for GO term: {prediction['go_term_id']}, sequence: {prediction['target_sequence_id']}"
                )

            # Confirmar los cambios
            self.session.commit()
            self.logger.info("Predictions successfully stored.")

        except Exception as e:
            self.logger.error(f"Error storing predictions: {e}")
            self.session.rollback()

