from itertools import combinations

from protein_metamorphisms_is.tasks.queue import QueueTaskInitializer
from protein_metamorphisms_is.sql.model import ProteinGOTermAssociations, GOTerm, GOPairs, GOResultsPairwise, SequenceEmbeddingGOAnnotationTransfer
from goatools import obo_parser
from goatools.semantic import TermCounts, resnik_sim, min_branch_length, get_info_content
from goatools.anno.gaf_reader import GafReader


class GoMultifunctionalityMetrics(QueueTaskInitializer):
    def __init__(self, conf):
        """
        Initializes the GoMetrics class with configuration settings for GO analysis.

        Args:
            conf (dict): Configuration dictionary containing paths to OBO and GAF files.
        """
        super().__init__(conf)
        self.go = obo_parser.GODag(conf['obo'])
        self.logger.info("GoMetrics instance created")
        self.reference_attribute = "go_terms_per_protein"

    def load_annotations_from_gaf(self, gaf_path):
        """Loads GO annotations from a GAF file."""
        return GafReader(gaf_path).read_gaf()

    def enqueue(self):
        """
        Loads GO terms for both proteins and clusters and publishes tasks for each GO term pair.
        """
        pairs = self.load_pairs()
        for pair in pairs:
            self.publish_task({'pair': pair})

    def load_go_terms_per_protein(self):
        """
        Carga y organiza los términos GO por proteína por categoría (P, C, F).

        Returns:
            dict: Diccionario que mapea los nombres de las proteínas a sus términos GO categorizados por tipo.
        """
        session = self.session
        go_terms_per_protein = {}

        try:
            # Realizar una consulta única para cargar todas las asociaciones y términos
            associations = session.query(ProteinGOTermAssociations, GOTerm).join(
                GOTerm, ProteinGOTermAssociations.go_id == GOTerm.go_id).all()

            for association, go_term in associations:
                if association.protein_entry_name not in go_terms_per_protein:
                    go_terms_per_protein[association.protein_entry_name] = {"P": [], "C": [], "F": []}

                category = go_term.category
                if category in go_terms_per_protein[association.protein_entry_name]:
                    go_terms_per_protein[association.protein_entry_name][category].append(go_term.go_id)
                else:
                    self.logger.warning(f"Categoría GO desconocida '{category}' para el término {go_term.go_id}.")
        except Exception as e:
            self.logger.error(f"Error cargando los términos GO por proteína: {e}")
            session.rollback()
        finally:
            session.close()

        return go_terms_per_protein

    def load_go_terms_per_cluster(self):
        """
        Carga y organiza los términos GO por clúster basado en la tabla de anotaciones transferidas `SequenceEmbeddingGOAnnotationTransfer`.

        Returns:
            dict: Diccionario que mapea IDs de clústeres a sus términos GO categorizados por tipo (P, C, F).
        """
        session = self.session
        go_terms_per_cluster = {}

        try:
            cluster_annotations = session.query(SequenceEmbeddingGOAnnotationTransfer, GOTerm).join(
                GOTerm, SequenceEmbeddingGOAnnotationTransfer.go_id == GOTerm.go_id).all()

            for annotation, go_term in cluster_annotations:
                cluster_id = annotation.source_cluster_id

                if cluster_id not in go_terms_per_cluster:
                    go_terms_per_cluster[cluster_id] = {"P": [], "C": [], "F": []}

                category = go_term.category
                if category in go_terms_per_cluster[cluster_id]:
                    go_terms_per_cluster[cluster_id][category].append(go_term.go_id)
                else:
                    self.logger.warning(f"Categoría GO desconocida '{category}' para el término {go_term.go_id}.")
        except Exception as e:
            self.logger.error(f"Error cargando los términos GO por clúster: {e}")
            session.rollback()
        finally:
            session.close()

        return go_terms_per_cluster

    def load_all_go_terms(self):
        """
        Combina los términos GO por proteína y por clúster en una única estructura de datos.

        Returns:
            dict: Diccionario que contiene los términos GO para proteínas y clústeres categorizados por tipo.
        """
        go_terms_per_protein = self.load_go_terms_per_protein()
        print(len(go_terms_per_protein))
        go_terms_per_cluster = self.load_go_terms_per_cluster()
        print(len(go_terms_per_cluster))

        # Unificar términos en una estructura común
        all_go_terms = {}

        # Agregar los términos de proteínas
        for protein, terms_by_category in go_terms_per_protein.items():
            all_go_terms[protein] = terms_by_category

        # Agregar los términos de clústeres
        for cluster_id, terms_by_category in go_terms_per_cluster.items():
            all_go_terms[f"Cluster_{cluster_id}"] = terms_by_category

        return all_go_terms

    def load_pairs(self):
        """
        Genera pares de términos GO para cada categoría (P, C, F) en proteínas y clústeres.

        Returns:
            list: Una lista de diccionarios que representan pares de términos GO por categoría y fuente.
        """
        all_go_terms = self.load_all_go_terms()
        all_pairs = []

        for source, terms_by_category in all_go_terms.items():
            for category, terms in terms_by_category.items():
                sorted_terms = sorted(terms)
                for go_term_1, go_term_2 in combinations(sorted_terms, 2):
                    all_pairs.append({
                        'go_term_1': go_term_1,
                        'go_term_2': go_term_2,
                        'source': source,
                        'category': category
                    })
                    self.logger.info(f"Generated pair: Source={source}, Category={category}, GO1={go_term_1}, GO2={go_term_2}")

        self.logger.info(f"Total pairs generated: {len(all_pairs)}")
        return all_pairs



    def process(self, data):
        """
        Processes the GO term pairs and calculates metrics for each pair.

        Args:
            data (dict): Dictionary containing 'pair' with the relevant GO term IDs and metadata.

        Returns:
            dict: A dictionary containing the calculated metrics for each GO term pair.
        """
        # Obtener el par de términos y la información del diccionario de datos
        pair = data['pair']
        go_term_1_id = pair['go_term_1']
        go_term_2_id = pair['go_term_2']

        # Calcular el contenido de información para cada término GO
        # term_counts = TermCounts(self.go, self.annotations)
        try:
            # ic1 = get_info_content(go_term_1_id, term_counts)
            # ic2 = get_info_content(go_term_2_id, term_counts)

            # Calcular la distancia de Resnik entre los términos GO
            # resnik = resnik_sim(go_term_1_id, go_term_2_id, self.go, term_counts)
            # Calcular la longitud mínima de rama entre los términos GO
            mbl = min_branch_length(go_term_1_id, go_term_2_id, self.go, branch_dist=None)
            # Guardar los resultados en un diccionario con las métricas calculadas

            result = {
                'go_term_1_id': go_term_1_id,
                'go_term_2_id': go_term_2_id,
                # 'information_content_1': ic1,
                # 'information_content_2': ic2,
                # 'resnik_distance': resnik,
                'minimum_branch_length': mbl
            }

            # Imprimir para verificar el cálculo
            self.logger.info(f"Calculated metrics for pair: {result}")
            return result

        except KeyError as e:
            self.logger.error(f"Error calculating metrics for terms {go_term_1_id} and {go_term_2_id}: {e}")
            return None

    def store_entry(self, record):
        """
        Stores or updates the calculated metrics into the database.

        Args:
            record (dict): Dictionary containing the calculated metrics for a specific pair of GO terms.
        """
        session = self.session
        go_term_1_id = record['go_term_1_id']
        go_term_2_id = record['go_term_2_id']

        try:
            # Buscar si el par ya existe en la tabla GOPairs
            existing_pair = session.query(GOPairs).filter(
                ((GOPairs.go_term_1_id == go_term_1_id) & (GOPairs.go_term_2_id == go_term_2_id)) |
                ((GOPairs.go_term_1_id == go_term_2_id) & (GOPairs.go_term_2_id == go_term_1_id))
            ).first()

            if existing_pair:
                # Si el par existe, usar el ID del par existente
                pair_id = existing_pair.id
                self.logger.info(f"Pair exists: {go_term_1_id}, {go_term_2_id} with Pair ID: {pair_id}")
            else:
                # Si el par no existe, crearlo
                new_pair = GOPairs(go_term_1_id=go_term_1_id, go_term_2_id=go_term_2_id)
                session.add(new_pair)
                session.flush()  # Asegura que new_pair.id esté disponible
                pair_id = new_pair.id
                self.logger.info(f"New pair created: {go_term_1_id}, {go_term_2_id} with Pair ID: {pair_id}")

            # Verificar si ya existen resultados para este par
            existing_result = session.query(GOResultsPairwise).filter_by(pair_id=pair_id).first()

            if existing_result:
                # Si ya existen resultados, actualizarlos
                # existing_result.information_content_1 = record['information_content_1']
                # existing_result.information_content_2 = record['information_content_2']
                # existing_result.resnik_distance = record['resnik_distance']
                existing_result.minimum_branch_length = record['minimum_branch_length']
                self.logger.info(f"Updated results for Pair ID: {pair_id}")
            else:
                # Crear nuevos resultados si no existen
                new_result = GOResultsPairwise(
                    pair_id=pair_id,
                    # information_content_1=record['information_content_1'],
                    # information_content_2=record['information_content_2'],
                    # resnik_distance=record['resnik_distance'],
                    minimum_branch_length=record['minimum_branch_length']
                )
                session.add(new_result)
                self.logger.info(f"Stored new results for Pair ID: {pair_id}")

            # Confirmar la transacción
            session.commit()

        except Exception as e:
            self.logger.error(f"Error saving GO term relationships: {e}")
            session.rollback()
