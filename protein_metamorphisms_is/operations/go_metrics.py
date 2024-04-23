from protein_metamorphisms_is.operations.base.operator import OperatorBase
from protein_metamorphisms_is.sql.model import ProteinGOTermAssociation, GOTerm, GOTermRelationship
from goatools import obo_parser
from goatools.semantic import TermCounts, resnik_sim, min_branch_length, get_info_content
from goatools.anno.gaf_reader import GafReader


class GoMetrics(OperatorBase):
    def __init__(self, conf):
        super().__init__(conf)
        self.go = obo_parser.GODag(conf['obo'])
        self.annotations = self.load_annotations_from_gaf(conf['go_annotation_file'])
        self.logger.info("GoMetrics instance created")

    def load_annotations_from_gaf(self, gaf_path):
        gaf = GafReader(gaf_path).read_gaf()
        return gaf

    def start(self):
        try:
            self.logger.info("Starting GO Metrics calculation")
            go_terms_per_protein = self.load_go_terms_per_protein()
            self.calculate_metrics(go_terms_per_protein)

            self.logger.info("GO Metrics calculation process completed successfully")

        except Exception as e:
            self.logger.error(f"Error during GO Metrics calculation process: {e}")
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
        term_counts = TermCounts(self.go, self.annotations)
        missing_go_terms = set()  # Para registrar los términos GO faltantes

        for protein, terms_per_category in go_terms_per_protein.items():
            for category, terms in terms_per_category.items():
                for i, go_term_1 in enumerate(terms):
                    for go_term_2 in terms[i + 1:]:
                        if go_term_1 != go_term_2:
                            # Asegurar orden alfabético para la consistencia
                            go_term_1_id, go_term_2_id = sorted([go_term_1, go_term_2])

                            # Verificar si la relación ya existe en la base de datos
                            exists = self.session.query(GOTermRelationship).filter(
                                ((GOTermRelationship.go_term_1_id == go_term_1_id) & (
                                        GOTermRelationship.go_term_2_id == go_term_2_id)) |
                                ((GOTermRelationship.go_term_1_id == go_term_2_id) & (
                                        GOTermRelationship.go_term_2_id == go_term_1_id))
                            ).first()

                            if exists:
                                continue

                            try:
                                ic1 = get_info_content(go_term_1_id, term_counts)
                                ic2 = get_info_content(go_term_2_id, term_counts)
                                resnik = resnik_sim(go_term_1_id, go_term_2_id, self.go, term_counts)
                                mbl = min_branch_length(go_term_1_id, go_term_2_id, self.go, branch_dist=None)
                            except KeyError as e:
                                missing_go_terms.add(str(e))
                                continue

                            self.save_go_term_relationship(go_term_1_id, go_term_2_id, ic1, ic2, resnik, mbl)

        if missing_go_terms:
            self.logger.error(f"Missing GO terms: {missing_go_terms}")

    def save_go_term_relationship(self, go_term_1_id, go_term_2_id, information_content_1, information_content_2,
                                  resnik_distance,
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
