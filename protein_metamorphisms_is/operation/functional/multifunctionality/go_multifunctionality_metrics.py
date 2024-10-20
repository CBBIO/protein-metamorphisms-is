from itertools import combinations

from sqlalchemy.orm import aliased

from protein_metamorphisms_is.tasks.queue import QueueTaskInitializer
from protein_metamorphisms_is.sql.model.model import GOTerm, GOPairs, \
    GOAnnotation, GOPairsTerms, GOResultsPairwise, ProteinGoPair
from goatools import obo_parser
from goatools.semantic import min_branch_length
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

    def load_filtered_annotations(self):
        """
        Carga las anotaciones de GO, filtrando aquellas que tienen un 'protein_entry'.

        Returns:
            list: Lista de anotaciones filtradas.
        """
        session = self.session
        try:
            # Filtrar solo las anotaciones con 'protein_entry'
            annotations = session.query(GOAnnotation).filter(GOAnnotation.protein_entry_name.isnot(None)).all()
            self.logger.info(f"Loaded {len(annotations)} annotations with 'protein_entry'.")
            return annotations
        except Exception as e:
            self.logger.error(f"Error loading filtered annotations: {e}")
            session.rollback()
        finally:
            session.close()

    def load_pairs(self):
        """
        Genera pares de términos GO para cada categoría (P, C, F) en proteínas.

        Returns:
            list: Una lista de diccionarios que representan pares de términos GO por categoría y proteína.
        """
        session = self.session
        all_pairs = []

        try:
            # Consultar las anotaciones con términos GO y proteínas
            annotations = (
                session.query(GOAnnotation, GOTerm)
                .join(GOTerm, GOAnnotation.go_id == GOTerm.go_id)
                .filter(GOAnnotation.protein_entry_name.isnot(None))
                .order_by(GOAnnotation.protein_entry_name, GOAnnotation.is_transferred.asc(),
                          GOAnnotation.distance.asc())
                .all()
            )

            # Organizar los términos GO por categoría y proteína
            go_terms_by_protein = {}
            for annotation, go_term in annotations:
                protein = annotation.protein_entry_name
                category = go_term.category

                # Inicializar la estructura para la proteína si no existe
                if protein not in go_terms_by_protein:
                    go_terms_by_protein[protein] = {'P': [], 'C': [], 'F': []}

                # Añadir el término GO a la categoría correspondiente
                if category in go_terms_by_protein[protein]:
                    go_terms_by_protein[protein][category].append(go_term.go_id)

            # Debug: Imprimir los términos GO organizados
            print(go_terms_by_protein)

            # Generar pares de términos GO por cada proteína y categoría
            for protein, categories in go_terms_by_protein.items():
                self.logger.info('cat', categories)
                for category, terms in categories.items():
                    sorted_terms = sorted(set(terms))
                    for go_term_1, go_term_2 in combinations(sorted_terms, 2):
                        all_pairs.append({
                            'go_term_1': go_term_1,
                            'go_term_2': go_term_2,
                            'category': category,
                            'protein': protein
                        })
                        self.logger.info(
                            f"Generated pair: Protein={protein}, Category={category}, GO1={go_term_1}, GO2={go_term_2}")


            self.logger.info(f"Total pairs generated with 'protein_entry': {len(all_pairs)}")
            return all_pairs

        except Exception as e:
            self.logger.error(f"Error loading GO pairs with 'protein_entry': {e}")
            session.rollback()

        finally:
            session.close()

    def process(self, data):
        """
        Processes the GO term pairs and calculates metrics for each pair.

        Args:
            data (dict): Dictionary containing 'pair' with the relevant GO term IDs and metadata.

        Returns:
            dict: A dictionary containing the calculated metrics for each GO term pair.
        """
        pair = data['pair']
        go_term_1_id = pair['go_term_1']
        go_term_2_id = pair['go_term_2']
        protein = pair['protein']

        # Calcular el contenido de información para cada término GO
        try:
            # term_counts = TermCounts(self.go, self.annotations)
            # ic1 = get_info_content(go_term_1_id, term_counts)
            # ic2 = get_info_content(go_term_2_id, term_counts)

            # Calcular la distancia de Resnik entre los términos GO
            # resnik = resnik_sim(go_term_1_id, go_term_2_id, self.go, term_counts)
            # Calcular la longitud mínima de rama entre los términos GO
            mbl = min_branch_length(go_term_1_id, go_term_2_id, self.go, branch_dist=None)

            # Guardar los resultados en un diccionario con las métricas calculadas
            result = {
                'protein': protein,
                'go_term_1_id': go_term_1_id,
                'go_term_2_id': go_term_2_id,
                # 'information_content_1': ic1,
                # 'information_content_2': ic2,
                # 'resnik_distance': resnik,
                'minimum_branch_length': mbl
            }

            self.logger.info(f"Calculated metrics for pair: {result}")
            return result

        except Exception as e:
            self.logger.error(f"Error calculating metrics for terms {go_term_1_id} and {go_term_2_id}: {e}")
            return None

    def store_entry(self, record):
        session = self.session
        try:
            # Extraer información del registro
            go_term_1_id = record['go_term_1_id']
            go_term_2_id = record['go_term_2_id']
            protein = record['protein']
            minimum_branch_length = record['minimum_branch_length']
            information_content_1 = record.get('information_content_1', None)
            information_content_2 = record.get('information_content_2', None)
            resnik_distance = record.get('resnik_distance', None)

            # Aliases para las dos instancias de GOPairsTerms
            gpt1 = aliased(GOPairsTerms)
            gpt2 = aliased(GOPairsTerms)

            # Paso 1: Verificar si el par de GO ya existe
            go_pair = session.query(GOPairs).join(gpt1, GOPairs.id == gpt1.id).join(
                gpt2, GOPairs.id == gpt2.id
            ).filter(
                gpt1.go_term_id == go_term_1_id,
                gpt2.go_term_id == go_term_2_id
            ).first()

            # Paso 2: Si el par no existe, crearlo
            if not go_pair:
                go_pair = GOPairs()
                session.add(go_pair)
                session.flush()  # Obtener el ID generado automáticamente

                # Almacenar los términos en GOPairsTerms
                go_pair_terms = [
                    GOPairsTerms(id=go_pair.id, go_term_id=go_term_1_id),
                    GOPairsTerms(id=go_pair.id, go_term_id=go_term_2_id)
                ]
                session.add_all(go_pair_terms)

            # Paso 3: Verificar si la relación en ProteinGoPair ya existe
            existing_protein_go_pair = session.query(ProteinGoPair).filter_by(
                go_pair_id=go_pair.id, protein_entry_name=protein
            ).first()

            # Si no existe, agregar la relación
            if not existing_protein_go_pair:
                protein_go_pair = ProteinGoPair(go_pair_id=go_pair.id, protein_entry_name=protein)
                session.add(protein_go_pair)

            # Paso 4: Almacenar los resultados en GOResultsPairwise
            results_entry = GOResultsPairwise(
                pair_id=go_pair.id,
                information_content_1=information_content_1,
                information_content_2=information_content_2,
                resnik_distance=resnik_distance,
                minimum_branch_length=minimum_branch_length
            )
            session.add(results_entry)

            # Commit de la transacción
            session.commit()
            self.logger.info(
                f"Stored entries for GO pair: {go_term_1_id} and {go_term_2_id} with pair ID {go_pair.id} and protein {protein}")

        except Exception as e:
            self.logger.error(f"Error storing entry for {go_term_1_id} and {go_term_2_id}: {e}")
            session.rollback()

        finally:
            session.close()
