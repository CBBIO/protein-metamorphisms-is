from itertools import combinations

from goatools.base import get_godag
from goatools.godag.go_tasks import get_go2ancestors

from protein_metamorphisms_is.sql.model.entities.go_annotation.go_annotation import ProteinGOTermAnnotation
from protein_metamorphisms_is.sql.model.entities.go_annotation.go_term import GOTerm
from protein_metamorphisms_is.sql.model.operational.functional.group import GOTermPair, GOTermPairEntry, \
    GOTermPairProtein
from protein_metamorphisms_is.sql.model.operational.functional.result import GOTermPairResult
from protein_metamorphisms_is.tasks.queue import QueueTaskInitializer

from goatools.anno.gaf_reader import GafReader


class GoMultifunctionalityMetrics(QueueTaskInitializer):
    def __init__(self, conf):
        """
        Initializes the GoMetrics class with configuration settings for GO analysis.

        Args:
            conf (dict): Configuration dictionary containing paths to OBO and GAF files.
        """
        super().__init__(conf)
        self.obo_path = self.conf.get('obo', '../data/go-basic.obo')
        self.go = get_godag(self.obo_path, optional_attrs='relationship')
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

    def load_pairs(self):
        """
        Generates unique pairs of GO terms per category (P, C, F) in proteins
        and associates each pair with a list of proteins where they appear.

        Returns:
            list: A list of dictionaries representing unique pairs of GO terms
                  per category, along with the associated proteins.
        """
        session = self.session
        pairs_dict = {}

        try:
            # Consultar las anotaciones con términos GO y proteínas
            annotations = (
                session.query(ProteinGOTermAnnotation, GOTerm)
                .join(GOTerm, ProteinGOTermAnnotation.go_id == GOTerm.go_id)
                .filter(ProteinGOTermAnnotation.id.isnot(None))
                .order_by(ProteinGOTermAnnotation.id)
                .all()
            )

            # Organizar los términos GO por categoría y proteína
            go_terms_by_protein = {}
            for annotation, go_term in annotations:
                protein = annotation.protein_id
                category = go_term.category

                # Inicializar la estructura para la proteína si no existe
                if protein not in go_terms_by_protein:
                    go_terms_by_protein[protein] = {'P': [], 'C': [], 'F': []}

                # Añadir el término GO a la categoría correspondiente
                if category in go_terms_by_protein[protein]:
                    go_terms_by_protein[protein][category].append(go_term.go_id)

            # Generar pares únicos de términos GO por categoría
            for protein, categories in go_terms_by_protein.items():
                for category, terms in categories.items():
                    sorted_terms = sorted(set(terms))
                    for go_term_1, go_term_2 in combinations(sorted_terms, 2):
                        pair_key = (go_term_1, go_term_2, category)
                        if pair_key not in pairs_dict:
                            pairs_dict[pair_key] = {'proteins': []}
                        pairs_dict[pair_key]['proteins'].append(protein)

            # Convertir el diccionario de pares en una lista estructurada
            all_pairs = [
                {
                    'go_term_1': go_term_1,
                    'go_term_2': go_term_2,
                    'category': category,
                    'proteins': sorted(set(data['proteins']))
                }
                for (go_term_1, go_term_2, category), data in pairs_dict.items()
            ]

            self.logger.info(f"Total unique pairs generated: {len(all_pairs)}")
            return all_pairs

        except Exception as e:
            self.logger.error(f"Error loading GO term pairs: {e}")
            return []

    def process(self, data):
        """
        Processes the GO term pairs and calculates metrics for each pair.

        Args:
            data (dict): Dictionary containing 'pair' with the relevant GO term IDs and metadata.

        Returns:
            dict: A dictionary containing the calculated metrics for each GO term pair.
        """

        # Calcular el contenido de información para cada término GO
        try:
            pair = data['pair']
            go_term_1_id = pair['go_term_1']
            go_term_2_id = pair['go_term_2']
            proteins = pair['proteins']

            # Calcular la MBL utilizando las relaciones opcionales
            mbl = calculate_mbl_with_relationships(go_term_1_id, go_term_2_id, self.go)

            # Guardar los resultados en un diccionario con las métricas calculadas
            result = {
                'proteins': proteins,
                'go_term_1_id': go_term_1_id,
                'go_term_2_id': go_term_2_id,
                'minimum_branch_length': mbl
            }

            self.logger.info(f"Calculated metrics for pair: {result}")
            return result

        except Exception as e:
            self.logger.error(f"Error calculating metrics for terms {go_term_1_id} and {go_term_2_id}: {e}")
            return None

    def store_entry(self, record):
        """
        Stores a GOTermPair, its entries, and the related proteins.

        Args:
            record (dict): Dictionary containing the GO term pair data and the associated proteins.
        """
        session = self.session
        try:
            # Extraer datos del registro
            go_term_1_id = record['go_term_1_id']
            go_term_2_id = record['go_term_2_id']
            proteins = record['proteins']
            mbl = record.get('minimum_branch_length', None)

            # Crear un nuevo GOTermPair
            go_term_pair = GOTermPair()
            session.add(go_term_pair)
            session.flush()  # Obtener el ID del GOTermPair

            # Crear entradas en GOTermPairEntry para los dos términos GO
            for go_term_id in [go_term_1_id, go_term_2_id]:
                term_entry = GOTermPairEntry(go_term_pair_id=go_term_pair.id, go_term_id=go_term_id)
                session.add(term_entry)

            # Asociar proteínas al GOTermPair mediante GOTermPairProtein
            for protein_id in proteins:
                association = session.query(GOTermPairProtein).filter_by(
                    go_term_pair_id=go_term_pair.id, protein_id=protein_id
                ).first()
                if not association:
                    association = GOTermPairProtein(go_term_pair_id=go_term_pair.id, protein_id=protein_id)
                    session.add(association)

            # Crear el resultado asociado al GOTermPair
            term_result = session.query(GOTermPairResult).filter_by(go_term_pair_id=go_term_pair.id).first()
            if not term_result:
                term_result = GOTermPairResult(go_term_pair_id=go_term_pair.id, mbl=mbl)
                session.add(term_result)

            # Confirmar los cambios
            session.commit()

        except Exception as e:
            self.logger.error(f"Error storing entry: {e}")
            session.rollback()
        finally:
            session.close()


def calculate_mbl_with_relationships(go_id1, go_id2, godag):
    # Crear el subgrafo con todas las relaciones
    go2ancestors = get_go2ancestors(set(godag.values()), relationships={"is_a", "part_of", "regulates"})

    # Obtener todos los ancestros de ambos términos
    ancestors1 = get_all_ancestors(go_id1, go2ancestors) | {go_id1}
    ancestors2 = get_all_ancestors(go_id2, go2ancestors) | {go_id2}

    # Encontrar los ancestros comunes
    common_ancestors = ancestors1.intersection(ancestors2)
    if not common_ancestors:
        print("No hay ancestros comunes entre los términos dados.")
        return None

    # Calcular la distancia mínima a los ancestros comunes
    min_distance = float("inf")
    for ancestor in common_ancestors:
        distance = abs(godag[go_id1].depth + godag[go_id2].depth - 2 * godag[ancestor].depth)
        min_distance = min(min_distance, distance)

    print(f"Minimum Branch Length (MBL) entre {go_id1} y {go_id2}: {min_distance}")
    return min_distance


def get_all_ancestors(go_id, go2ancestors):
    """Devuelve todos los ancestros del término GO incluyendo relaciones opcionales."""
    return go2ancestors.get(go_id, set())
