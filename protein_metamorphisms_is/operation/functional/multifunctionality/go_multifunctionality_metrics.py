from itertools import combinations

from protein_metamorphisms_is.tasks.queue import QueueTaskInitializer
from protein_metamorphisms_is.sql.model import ProteinGOTermAssociations, GOTerm, GOPairs, GOResultsPairwise
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
        self.annotations = self.load_annotations_from_gaf(conf['go_annotation_file'])
        self.logger.info("GoMetrics instance created")
        self.reference_attribute = "go_terms_per_protein"

    def load_annotations_from_gaf(self, gaf_path):
        """Loads GO annotations from a GAF file."""
        return GafReader(gaf_path).read_gaf()

    def enqueue(self):
        """
        Loads GO terms per protein and publishes tasks for each protein and its GO term categories.
        """
        go_terms_per_protein = self.load_go_terms_per_protein()
        for protein, terms_per_category in go_terms_per_protein.items():
            self.publish_task({'protein': protein, 'terms_per_category': terms_per_category})

    def load_go_terms_per_protein(self):
        """
        Loads and organizes GO terms per protein by category (P, C, F).

        Returns:
            dict: Dictionary mapping protein entry names to their GO terms categorized by type.
        """
        session = self.session
        go_terms_per_protein = {}

        try:
            # Perform a single query to load all associations and terms
            associations = session.query(ProteinGOTermAssociations, GOTerm).join(
                GOTerm, ProteinGOTermAssociations.go_id == GOTerm.go_id).all()

            for association, go_term in associations:
                if association.protein_entry_name not in go_terms_per_protein:
                    go_terms_per_protein[association.protein_entry_name] = {"P": [], "C": [], "F": []}

                category = go_term.category
                if category in go_terms_per_protein[association.protein_entry_name]:
                    go_terms_per_protein[association.protein_entry_name][category].append(go_term.go_id)
                else:
                    self.logger.warning(f"Unknown GO category '{category}' for term {go_term.go_id}.")
        except Exception as e:
            self.logger.error(f"Error loading GO terms per protein: {e}")
            session.rollback()
        finally:
            session.close()

        print(go_terms_per_protein)
        return go_terms_per_protein

    def process(self, data):
        """
        Processes the GO term data for a given protein and computes metrics.

        Args:
            data (dict): Dictionary containing 'protein' and 'terms_per_category'.

        Returns:
            list: A list of metrics for each GO term pair.
        """
        term_counts = TermCounts(self.go, self.annotations)
        protein = data['protein']
        terms_per_category = data['terms_per_category']
        metrics = []

        for category, terms in terms_per_category.items():
            self.logger.debug(f'Processing category {category} for protein {protein} with terms {terms}')
            metrics.extend(self.calculate_pairwise_metrics(terms, term_counts, protein, category))

        return metrics

    def calculate_pairwise_metrics(self, terms, term_counts, protein, category):
        """
        Calcula métricas de pares de términos GO utilizando combinatoria sin repetición.

        Args:
            terms (list): Lista de términos GO de una categoría específica.
            term_counts (TermCounts): Conteos de términos precomputados para los términos GO.
            protein (str): Nombre del protein.
            category (str): Categoría de los términos GO (P, C, F).

        Returns:
            list: Lista de diccionarios que contienen las métricas calculadas.
        """
        metrics = []

        # Generar combinaciones únicas de los términos sin repeticiones
        for go_term_1, go_term_2 in combinations(terms, 2):
            # Aplicar la regla de orden: go_term_1 debe ser mayor que go_term_2
            if go_term_1 > go_term_2:
                # Hacer el cálculo si no existe la relación
                if not self.pairwise_relationship_exists(go_term_1, go_term_2):
                    try:
                        ic1 = get_info_content(go_term_1, term_counts)
                        ic2 = get_info_content(go_term_2, term_counts)
                        resnik = resnik_sim(go_term_1, go_term_2, self.go, term_counts)
                        mbl = min_branch_length(go_term_1, go_term_2, self.go, branch_dist=None)

                        # Guardar las métricas calculadas
                        metrics.append({
                            'protein_entry_name': protein,
                            'category': category,
                            'go_term_1_id': go_term_1,
                            'go_term_2_id': go_term_2,
                            'information_content_1': ic1,
                            'information_content_2': ic2,
                            'resnik_distance': resnik,
                            'minimum_branch_length': mbl
                        })
                    except KeyError as e:
                        self.logger.error(f"Missing GO terms: {str(e)}")
        return metrics

    def pairwise_relationship_exists(self, go_term_1_id, go_term_2_id):
        """
        Checks if a pairwise relationship between two GO terms already exists in the database.

        Args:
            go_term_1_id (str): GO term ID of the first term.
            go_term_2_id (str): GO term ID of the second term.

        Returns:
            bool: True if relationship exists, False otherwise.
        """
        return self.session.query(GOPairs).filter(
            ((GOPairs.go_term_1_id == go_term_1_id) & (GOPairs.go_term_2_id == go_term_2_id)) |
            ((GOPairs.go_term_1_id == go_term_2_id) & (GOPairs.go_term_2_id == go_term_1_id))
        ).first() is not None

    def store_entry(self, record):
        """
        Stores the calculated metrics into the database.

        Args:
            record (list): List of metrics for a specific protein and GO term pairs.
        """
        session = self.session

        try:
            for entry in record:
                # Verificar si el par ya existe en la tabla
                existing_pair = session.query(GOPairs).filter_by(
                    go_term_1_id=entry['go_term_1_id'],
                    go_term_2_id=entry['go_term_2_id']
                ).first()

                if existing_pair:
                    new_pair_id = existing_pair.id
                else:
                    # Crear un nuevo par si no existe
                    new_pair = GOPairs(
                        go_term_1_id=entry['go_term_1_id'],
                        go_term_2_id=entry['go_term_2_id']
                    )
                    session.add(new_pair)
                    session.flush()  # Asegura que new_pair.id esté disponible
                    new_pair_id = new_pair.id
                    if new_pair_id is None:
                        self.logger.error(
                            f"Failed to generate pair_id for terms: {entry['go_term_1_id']} and {entry['go_term_2_id']}")
                        continue

                # Crear la entrada del resultado con el ID del par
                result_entry = GOResultsPairwise(
                    pair_id=new_pair_id,
                    information_content_1=entry['information_content_1'],
                    information_content_2=entry['information_content_2'],
                    resnik_distance=entry['resnik_distance'],
                    minimum_branch_length=entry['minimum_branch_length']
                )

                session.add(result_entry)

            session.commit()
        except Exception as e:
            self.logger.error(f"Error saving GO term relationships: {e}")
            session.rollback()
