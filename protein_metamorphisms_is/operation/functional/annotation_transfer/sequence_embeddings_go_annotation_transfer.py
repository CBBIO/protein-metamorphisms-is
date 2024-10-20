from sqlalchemy import text

from protein_metamorphisms_is.sql.model.model import ClusterEntry, SequenceEmbedding, \
    GOAnnotation
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
            for cluster_entry in clusters:
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
                WITH target_cluster AS (
                    SELECT se.embedding
                    FROM clusters c
                    JOIN cluster_entries ce ON ce.cluster_id = c.id
                    JOIN sequences s ON s.id = ce.sequence_id
                    JOIN sequence_embeddings se ON se.sequence_id = s.id
                    WHERE c.id =  :cluster_id AND ce.is_representative = TRUE AND se.embedding_type_id = :embedding_type_id
                ),
                all_representatives AS (
                    SELECT 
                        c.id AS cluster_id,
                        se.embedding
                    FROM clusters c
                    JOIN cluster_entries ce ON ce.cluster_id = c.id
                    JOIN sequences s ON s.id = ce.sequence_id
                    JOIN sequence_embeddings se ON se.sequence_id = s.id
                    WHERE se.embedding_type_id = :embedding_type_id AND ce.is_representative = TRUE
                ),
                top_k_clusters AS (
                    SELECT 
                        ar.cluster_id,
                        l2_distance(tc.embedding, ar.embedding) AS distance
                    FROM 
                        target_cluster tc, 
                        all_representatives ar
                    ORDER BY 
                        distance ASC
                    LIMIT :k
                ),
                go_term_data AS (
                    SELECT 
                        ce.cluster_id,
                        ce.id AS cluster_entry_id,
                        tk.distance,
                        go_assoc.go_id AS go_term_id,
                        p.entry_name AS protein_entry,
                        ROW_NUMBER() OVER (PARTITION BY go_assoc.go_id ORDER BY tk.distance ASC) AS rn
                    FROM 
                        top_k_clusters tk
                    LEFT JOIN cluster_entries ce ON ce.cluster_id = tk.cluster_id
                    LEFT JOIN sequences s ON s.id = ce.sequence_id                
                    LEFT JOIN pdb_chains pc ON pc.sequence_id = ce.sequence_id    
                    LEFT JOIN pdb_references pr ON pr.id = pc.pdb_reference_id   
                    JOIN proteins p ON p.entry_name = pr.protein_entry_name
                    LEFT JOIN go_annotation go_assoc ON go_assoc.protein_entry_name = p.entry_name
                    LEFT JOIN go_terms go ON go.go_id = go_assoc.go_id
                    WHERE go_assoc.is_transferred = FALSE
                )
                SELECT 
                    cluster_id,
                    cluster_entry_id,
                    distance,
                    go_term_id,
                    protein_entry
                FROM 
                    go_term_data
                WHERE 
                    rn = 1  -- Seleccionar solo la fila con menor distancia por cada go_term_id
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
        Almacenar las entradas de predicción en la base de datos, verificando si ya existen.
        """
        try:
            self.logger.info("Storing predictions...")

            for prediction in predictions:
                # Crear una nueva instancia de GOAnnotation
                transfer_entry = GOAnnotation(
                    go_id=prediction['go_term_id'],
                    protein_entry_name=prediction['protein_entry'],
                    source_cluster_id=prediction['source_cluster_id'],
                    target_cluster_id=prediction['target_cluster_id'],
                    distance=prediction['distance'],
                    is_transferred=True,
                    embedding_type_id=prediction['embedding_type_id']
                )

                # Añadir a la sesión
                self.session.add(transfer_entry)
                self.logger.info(
                    f"Stored new entry for protein: {prediction['protein_entry']}, GO term: {prediction['go_term_id']}")

            # Confirmar transacción
            self.session.commit()
            self.logger.info("Predictions successfully stored.")

        except Exception as e:
            self.logger.error(f"Error storing predictions: {e}")
            self.session.rollback()

