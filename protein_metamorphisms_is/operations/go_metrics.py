from protein_metamorphisms_is.operations.base.operator import OperatorBase
from protein_metamorphisms_is.sql.model import PDBChains, Cluster, ProteinGOTermAssociation, GOTerm, GOTermRelationship
from pycdhit import cd_hit, read_clstr
from goatools import obo_parser
from goatools.semantic import TermCounts, resnik_sim, min_branch_length, deepest_common_ancestor, get_info_content

class GoMetrics(OperatorBase):
    def __init__(self, conf):
        super().__init__(conf)
        self.go = obo_parser.GODag('../data/go-basic.obo')
        self.logger.info("GoMetrics instance created")

    def start(self):
        try:
            self.logger.info("Starting CD-HIT clustering process")
            go_terms_per_protein = self.load_go_terms_per_protein()
            self.calculate_metrics(go_terms_per_protein)
            # self.explore_representatives()

            self.logger.info("Clustering process completed successfully")

        except Exception as e:
            self.logger.error(f"Error during clustering process: {e}")
            raise

    def load_go_terms_per_protein(self):
        session = self.session
        go_terms_per_protein = {}
        try:
            # Realizar una consulta que incluya la información de categoría de los términos GO
            # Asumiendo que tienes un modelo GOTerm que incluye un campo 'category'
            associations = session.query(ProteinGOTermAssociation, GOTerm).join(GOTerm,
                                                                                ProteinGOTermAssociation.go_id == GOTerm.go_id).all()
            for association, go_term in associations:
                if association.protein_entry_name not in go_terms_per_protein:
                    go_terms_per_protein[association.protein_entry_name] = {"P": [],
                                                                            "C": [],
                                                                            "F": []}

                # Categoriza el término GO basado en su categoría
                category = go_term.category
                if category in go_terms_per_protein[association.protein_entry_name]:
                    go_terms_per_protein[association.protein_entry_name][category].append(go_term.go_id)
                else:
                    self.logger.warning(f"Unknown GO category '{category}' for term {go_term.go_id}.")
        except Exception as e:
            self.logger.error(f"Error loading GO terms per protein: {e}")
            session.rollback()
            raise
        finally:
            session.close()

        return go_terms_per_protein

    def calculate_metrics(self, go_terms_per_protein):
        # Construir el diccionario completo de anotaciones de GO para todas las proteínas
        all_go_terms = {}
        for protein, terms_per_category in go_terms_per_protein.items():
            all_terms = []
            for category, terms in terms_per_category.items():
                all_terms.extend(terms)
            all_go_terms[protein] = all_terms
        term_counts = TermCounts(self.go, all_go_terms)

        # Inicializar TermCounts una sola vez con todas las anotaciones


        for protein, terms_per_category in go_terms_per_protein.items():
            for category, terms in terms_per_category.items():
                # Iterar sobre todos los pares únicos de términos GO como antes
                for i, go_term_1 in enumerate(terms):
                    for go_term_2 in terms[i + 1:]:
                        if go_term_1 != go_term_2:
                            # Calcular el ancestro común más profundo (DCA)
                            ic1 = get_info_content(go_term_1, term_counts)
                            ic2 = get_info_content(go_term_2, term_counts)
                            resnik = resnik_sim(go_term_1, go_term_2, self.go, term_counts)
                            # Calcular la longitud mínima de la rama (MBL)
                            mbl = min_branch_length(go_term_1, go_term_2, self.go, branch_dist=None)

                            self.save_go_term_relationship(go_term_1, go_term_2, ic1,ic2, resnik, mbl)

    def save_go_term_relationship(self, go_term_1_id, go_term_2_id, information_content_1,information_content_2, resnik_distance,
                                  minimum_branch_length):
        session = self.session
        try:
            new_relationship = GOTermRelationship(
                go_term_1_id=go_term_1_id,
                go_term_2_id=go_term_2_id,
                information_content_1=information_content_1,
                information_content_2=information_content_2,
                resnik_distance=resnik_distance,
                minimum_branch_length=minimum_branch_length
            )
            session.add(new_relationship)
            session.commit()
        except Exception as e:
            self.logger.error(f"Error saving GO term relationship: {e}")
            session.rollback()
