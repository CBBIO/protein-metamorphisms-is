import multiprocessing

from sqlalchemy import text
from protein_metamorphisms_is.tasks.queue import QueueTaskInitializer
from protein_metamorphisms_is.sql.model.model import ClusterGOMultifunctionalityMetrics
from goatools.semantic import min_branch_length


class GoAnnotationTransferAnalysis(QueueTaskInitializer):
    def __init__(self, conf):
        """
        Constructor de la clase que inicializa la configuración y el logger.
        """
        super().__init__(conf)
        self.logger.info("GoAnnotationTransferAnalysis initialized")

    def enqueue(self):
        """
        Enqueue a list of GO terms by category associated with each cluster.
        """
        try:
            # Consulta para obtener todos los términos GO agrupados por clúster y categoría
            # Cambiar SequenceEmbeddingGOAnnotationTransfer por ProteinGOAnnotation
            query = text("""
            SELECT
                t.source_cluster_id AS cluster_id,
                go.category AS category,
                ARRAY_AGG(t.go_id) AS go_terms
            FROM protein_go_annotation t
            JOIN go_terms go ON go.go_id = t.go_id
            WHERE t.is_transferred = TRUE
            GROUP BY t.source_cluster_id, go.category
            """)

            with self.engine.connect() as connection:
                go_cluster_data = connection.execute(query).fetchall()

            # Encolar una tarea para cada clúster con sus términos GO por categoría
            for row in go_cluster_data:
                task_data = {
                    'cluster_id': row['cluster_id'],
                    'category': row['category'],
                    'go_terms': row['go_terms']
                }
                self.logger.info(f"Enqueued task for cluster {row['cluster_id']} with category {row['category']}")
                with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
                    pool.map(self.publish_task, task_data)

        except Exception as e:
            self.logger.error(f"Error enqueuing tasks: {e}")

    def process(self, task_data):
        """
        Calculate min_branch_length (MBL) for each pair of GO terms in the cluster.
        """
        try:
            cluster_id = task_data['cluster_id']
            category = task_data['category']
            go_terms = task_data['go_terms']

            # Verificar que haya suficientes términos GO para calcular MBL
            if len(go_terms) < 2:
                self.logger.warning(f"Not enough GO terms to calculate MBL for cluster {cluster_id}, category {category}")
                return None

            # Cargar el grafo GO desde el archivo de configuración
            go_dag = self.load_go_dag()  # Asegúrate de que esto cargue el grafo GO correctamente

            # Calcular MBL para cada par de términos GO
            total_mbl = 0
            count = 0

            for i in range(len(go_terms)):
                for j in range(i + 1, len(go_terms)):
                    term1 = go_terms[i]
                    term2 = go_terms[j]

                    try:
                        mbl = min_branch_length(term1, term2, go_dag)
                        total_mbl += mbl
                        count += 1
                        self.logger.info(f"MBL between {term1} and {term2} in cluster {cluster_id}: {mbl}")

                    except KeyError as e:
                        self.logger.warning(f"Error calculating MBL between {term1} and {term2}: {e}")

            # Calcular la media de MBL como métrica de multifuncionalidad
            multifunctionality_score = total_mbl / count if count > 0 else 0

            # Devolver el resultado
            result = {
                'cluster_id': cluster_id,
                'category': category,
                'multifunctionality_score': multifunctionality_score
            }

            return result

        except Exception as e:
            self.logger.error(f"Error processing GO terms for cluster {task_data['cluster_id']}: {e}")
            return None

    def store_entry(self, record):
        """
        Store the calculated multifunctionality metric in the database.
        """
        try:
            cluster_id = record['cluster_id']
            category = record['category']
            multifunctionality_score = record['multifunctionality_score']

            # Crear un nuevo registro en la tabla `ClusterGOMultifunctionalityMetrics`
            new_metric = ClusterGOMultifunctionalityMetrics(
                cluster_id=cluster_id,
                category=category,
                multifunctionality_score=multifunctionality_score
            )

            # Añadir y confirmar la entrada en la sesión
            self.session.add(new_metric)
            self.session.commit()
            self.logger.info(f"Stored multifunctionality score for cluster {cluster_id}, category {category}.")

        except Exception as e:
            self.logger.error(f"Error storing multifunctionality metrics: {e}")
            self.session.rollback()
