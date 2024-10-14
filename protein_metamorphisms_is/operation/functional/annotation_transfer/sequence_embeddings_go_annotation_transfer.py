from sqlalchemy import text

from protein_metamorphisms_is.sql.model import ClusterEntry, SequenceEmbedding, SequenceEmbeddingGOAnnotationTransfer
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
                    self.logger.warning(
                        f"No embedding found for sequence ID {sequence.id} in cluster {cluster_entry.cluster_id}")

        except Exception as e:
            self.logger.error(f"Error enqueuing tasks: {e}")

    def process(self, task_data):
        try:
            cluster_id = task_data['cluster_id']

            # Definir la consulta con la nueva lógica y estructura
            query = text("""
            WITH representative_cluster AS (
                SELECT se.embedding
                FROM clusters c
                JOIN cluster_entries ce ON ce.cluster_id = c.id
                JOIN sequences s ON s.id = ce.sequence_id
                JOIN sequence_embeddings se ON se.sequence_id = s.id
                WHERE c.id = :cluster_id AND ce.is_representative = TRUE
            ),
            all_representatives AS (
                SELECT 
                    c.id AS cluster_id,
                    se.embedding
                FROM clusters c
                JOIN cluster_entries ce ON ce.cluster_id = c.id
                JOIN sequences s ON s.id = ce.sequence_id
                JOIN sequence_embeddings se ON se.sequence_id = s.id
                WHERE ce.is_representative = TRUE
            ),
            top_k_clusters AS (
                SELECT 
                    ar.cluster_id,
                    l2_distance(rc.embedding, ar.embedding) AS distance
                FROM 
                    representative_cluster rc, 
                    all_representatives ar
                ORDER BY 
                    distance ASC
                LIMIT :k
            ),
            go_term_data AS (
                SELECT 
                    ce.cluster_id,
                    ce.id AS cluster_entry_id,
                    ce.sequence_id,
                    s.sequence,                   
                    ce.is_representative,
                    ce.sequence_length,
                    ce.identity,
                    ce.created_at,
                    tk.distance,
                    pc.chains AS chain_id,        
                    pr.pdb_id AS pdb_reference,   
                    pr.method AS pdb_method,      
                    pr.resolution,               
                    p.entry_name AS protein_name,
                    p.organism,                  
                    p.gene_name,                 
                    p.description,               
                    go.go_id AS go_term_id,      
                    go.category AS go_category,  
                    go.description AS go_description,
                    ROW_NUMBER() OVER (PARTITION BY go.go_id ORDER BY tk.distance ASC) AS rn
                FROM 
                    top_k_clusters tk
                JOIN cluster_entries ce ON ce.cluster_id = tk.cluster_id
                JOIN sequences s ON s.id = ce.sequence_id                
                LEFT JOIN pdb_chains pc ON pc.sequence_id = ce.sequence_id    
                LEFT JOIN pdb_references pr ON pr.id = pc.pdb_reference_id   
                JOIN proteins p ON p.sequence_id = ce.sequence_id
                JOIN protein_go_term_association go_assoc ON go_assoc.protein_entry_name = p.entry_name 
                JOIN go_terms go ON go.go_id = go_assoc.go_id   
            )
            SELECT 
                distance,
                cluster_id,
                cluster_entry_id,
                sequence_id,
                sequence,
                is_representative,
                sequence_length,
                identity,
                created_at,
                chain_id,
                pdb_reference,
                pdb_method,
                resolution,
                protein_name,
                organism,
                gene_name,
                description,
                go_term_id,
                go_category,
                go_description
            FROM 
                go_term_data
            WHERE 
                rn = 1  
            ORDER BY 
                cluster_id, go_term_id;
            """)

            # Obtener una conexión SQLAlchemy usando el motor
            with self.engine.connect() as connection:
                # Ejecutar la consulta con el parámetro cluster_id
                try:
                    results = connection.execute(query, {'cluster_id': cluster_id, 'k': self.k}).fetchall()
                    self.logger.info(f"Query executed successfully for cluster_id: {cluster_id}")
                except Exception as e:
                    self.logger.error(f"Error executing query: {e}")
                    return []

            # Convertir los resultados a una lista de diccionarios
            predictions = []
            try:
                for result in results:
                    prediction = {
                        'go_term_id': result.go_term_id,
                        'source_cluster_id': cluster_id,
                        'target_cluster_id': result.cluster_id,
                        'distance': result.distance
                    }
                    predictions.append(prediction)
            except Exception as e:
                self.logger.error(f"Error processing results: {e}")
                return []

            return predictions  # Devolver la lista de predicciones para su almacenamiento

        except Exception as e:
            self.logger.error(f"Unexpected error in process: {e}")
            return []

    def store_entry(self, predictions):
        """
        Store the prediction entries in the database.
        """
        try:
            self.logger.info("Storing predictions...")

            # Recorrer los resultados y preparar la inserción
            for prediction in predictions:
                # Crear una nueva instancia de SequenceEmbeddingGOAnnotationTransfer con `distance`
                transfer_entry = SequenceEmbeddingGOAnnotationTransfer(
                    go_id=prediction['go_term_id'],
                    source_cluster_id=prediction['source_cluster_id'],
                    target_cluster_id=prediction['target_cluster_id'],
                    distance=prediction['distance']  # Almacenar la distancia en la BD
                )

                # Añadir la nueva entrada a la sesión y confirmar los cambios
                self.session.add(transfer_entry)

            # Confirmar la transacción
            self.session.commit()
            self.logger.info("Predictions successfully stored.")

        except Exception as e:
            self.logger.error(f"Error storing predictions: {e}")
            self.session.rollback()
