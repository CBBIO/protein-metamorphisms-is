import json
import os

import pandas as pd
from sqlalchemy import text
import h5py

from protein_metamorphisms_is.tasks.queue import QueueTaskInitializer


class EmbeddingLookUp(QueueTaskInitializer):
    def __init__(self, conf):
        """
        Constructor de la clase que inicializa la configuración y el logger.
        """
        super().__init__(conf)
        self.logger.info("EmbeddingLookUp initialized")
        self.h5_path = conf.get('fantasia_output_h5', 'embeddings_results.h5')
        self.max_distance = conf.get('max_distance', 3)  # Distancia máxima permitida para el cálculo

    def enqueue(self):
        """
        Enqueue a task for each embedding from HDF5.
        """
        try:
            self.logger.info(f"Reading embeddings from HDF5: {self.h5_path}")

            tasks = []
            with h5py.File(self.h5_path, "r") as h5file:
                for accession, accession_group in h5file.items():
                    for embedding_type, type_group in accession_group.items():
                        # Verificar si el dataset 'embedding' existe
                        if 'embedding' in type_group:
                            embedding = type_group['embedding'][:]  # Leer el dataset como array
                            embedding_type_id = int(embedding_type.split('_')[1])  # Extraer el ID del tipo

                            task_data = {
                                'accession': accession,
                                'embedding': embedding,
                                'embedding_type_id': embedding_type_id
                            }
                            tasks.append(task_data)
                            self.logger.info(
                                f"Enqueued task for accession {accession} and embedding type {embedding_type_id}")
                        else:
                            self.logger.warning(f"Dataset 'embedding' not found in {embedding_type}")

            for task in tasks:
                self.publish_task(task)

            self.logger.info(f"Enqueued {len(tasks)} tasks based on HDF5 embeddings.")

        except Exception as e:
            self.logger.error(f"Error enqueuing tasks from HDF5: {e}")
            raise

    def process(self, task_data):
        """
        Procesa cada tarea para calcular términos GO basados en la distancia de embeddings.
        """
        try:
            accession = task_data['accession'].removeprefix('accession_')
            embedding_type_id = int(task_data['embedding_type_id'])  # Convertir a tipo nativo de Python
            embedding = task_data['embedding']

            # Convertir el vector NumPy a una lista de Python y formatear como cadena para pgvector
            vector_string = "[" + ",".join(f"{float(v):.8f}" for v in embedding) + "]"

            # Construir la consulta SQL con el vector embebido
            query = text(f"""
            WITH target_embedding AS (
                SELECT :vector_string ::vector AS embedding
            )
            SELECT 
                s.sequence,
                (se.embedding <-> te.embedding) AS distance,
                p.id AS protein_id,
                p.gene_name AS gene_name,
                p.organism AS organism,
                pgo.go_id AS go_term_id,
                gt.description AS go_term_description
            FROM 
                sequence_embeddings se
                JOIN target_embedding te ON TRUE
                JOIN sequence s ON se.sequence_id = s.id
                JOIN protein p ON s.id = p.sequence_id
                LEFT JOIN protein_go_term_annotation pgo ON p.id = pgo.protein_id
                LEFT JOIN go_terms gt ON pgo.go_id = gt.go_id
            WHERE 
                se.embedding_type_id = :embedding_type_id
                AND (se.embedding <-> te.embedding) < :max_distance
            ORDER BY 
                distance ASC;
            """)

            self.logger.info(f"Executing query for accession {accession} and embedding type {embedding_type_id}.")
            # Ejecutar la consulta
            with self.engine.connect() as connection:
                results = connection.execute(query, {
                    'vector_string': vector_string,
                    'embedding_type_id': embedding_type_id,
                    'max_distance': float(self.max_distance)
                }).fetchall()

            if not results:
                self.logger.info(f"No results found for accession {accession}.")
                return []

            go_terms = []
            for row in results:
                go_terms.append({
                    'accession': accession,
                    'go_id': row.go_term_id,
                    'go_description': row.go_term_description,
                    'distance': row.distance,
                    'embedding_type_id': embedding_type_id,
                    'protein_id': row.protein_id,
                    'organism': row.organism,
                    'sequence': row.sequence
                })

            self.logger.info(
                f"Found {len(go_terms)} GO terms for accession {accession} and embedding type {embedding_type_id}.")
            return go_terms

        except Exception as e:
            self.logger.error(
                f"Error processing task for accession {accession} and embedding type {embedding_type_id}: {e}")
            raise

    def store_entry(self, go_terms):
        """
        Guarda los términos GO en un archivo CSV utilizando pandas.
        Filtra las entradas donde go_terms sea None o esté vacío.
        """
        if not go_terms:  # Verifica si go_terms es None o está vacío
            self.logger.info("No valid GO terms to store.")
            return

        try:
            output_csv_path = self.conf.get('fantasia_output_csv', 'results.csv')

            # Convertir go_terms a un DataFrame
            df = pd.DataFrame(go_terms)

            # Verificar y crear el directorio si no existe
            output_dir = os.path.dirname(output_csv_path)
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
                self.logger.info(f"Created directory: {output_dir}")

            # Guardar en CSV (agregar al archivo existente si ya existe)
            if not os.path.exists(output_csv_path):
                df.to_csv(output_csv_path, index=False)
                self.logger.info(f"Results successfully stored in {output_csv_path}.")
            else:
                df.to_csv(output_csv_path, mode='a', index=False, header=False)
                self.logger.info(f"Appended results to {output_csv_path}.")

        except Exception as e:
            self.logger.error(f"Error storing results in CSV: {e}")
            raise
