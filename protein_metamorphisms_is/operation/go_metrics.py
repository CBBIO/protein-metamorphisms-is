from protein_metamorphisms_is.tasks.queue import QueueTaskInitializer
from protein_metamorphisms_is.sql.model import ProteinGOTermAssociation, GOTerm, GOTermRelationship
from goatools import obo_parser
from goatools.semantic import TermCounts, resnik_sim, min_branch_length, get_info_content
from goatools.anno.gaf_reader import GafReader


class GoMetrics(QueueTaskInitializer):
    def __init__(self, conf):
        super().__init__(conf)
        self.go = obo_parser.GODag(conf['obo'])
        self.annotations = self.load_annotations_from_gaf(conf['go_annotation_file'])
        self.logger.info("GoMetrics instance created")
        self.reference_attribute = "go_terms_per_protein"

    def load_annotations_from_gaf(self, gaf_path):
        return GafReader(gaf_path).read_gaf()

    def enqueue(self):
        go_terms_per_protein = self.load_go_terms_per_protein()
        for protein, terms_per_category in go_terms_per_protein.items():
            self.publish_task({'protein': protein, 'terms_per_category': terms_per_category})

    def load_go_terms_per_protein(self):
        session = self.session
        go_terms_per_protein = {}
        try:
            associations = session.query(ProteinGOTermAssociation, GOTerm).join(GOTerm,
                                                                                ProteinGOTermAssociation.go_id == GOTerm.go_id).all()

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
            raise
        finally:
            session.close()
        return go_terms_per_protein

    def process(self, data):
        term_counts = TermCounts(self.go, self.annotations)
        protein = data['protein']
        terms_per_category = data['terms_per_category']
        metrics = []

        for category, terms in terms_per_category.items():
            self.logger.debug(f'Processing category {category} for protein {protein} with terms {terms}')
            metrics.extend(self.calculate_pairwise_metrics(terms, term_counts))

        return metrics

    def calculate_pairwise_metrics(self, terms, term_counts):
        metrics = []
        for i, go_term_1 in enumerate(terms):
            for go_term_2 in terms[i + 1:]:
                if go_term_1 != go_term_2:
                    go_term_1_id, go_term_2_id = sorted([go_term_1, go_term_2])
                    if not self.relationship_exists(go_term_1_id, go_term_2_id):
                        try:
                            ic1 = get_info_content(go_term_1_id, term_counts)
                            ic2 = get_info_content(go_term_2_id, term_counts)
                            resnik = resnik_sim(go_term_1_id, go_term_2_id, self.go, term_counts)
                            mbl = min_branch_length(go_term_1_id, go_term_2_id, self.go, branch_dist=None)
                            metrics.append({
                                'go_term_1_id': go_term_1_id,
                                'go_term_2_id': go_term_2_id,
                                'information_content_1': ic1,
                                'information_content_2': ic2,
                                'resnik_distance': resnik,
                                'minimum_branch_length': mbl
                            })
                        except KeyError as e:
                            self.logger.error(f"Missing GO terms: {str(e)}")
        return metrics

    def relationship_exists(self, go_term_1_id, go_term_2_id):
        exists = self.session.query(GOTermRelationship).filter(
            ((GOTermRelationship.go_term_1_id == go_term_1_id) & (GOTermRelationship.go_term_2_id == go_term_2_id)) |
            ((GOTermRelationship.go_term_1_id == go_term_2_id) & (GOTermRelationship.go_term_2_id == go_term_1_id))
        ).first()
        return exists is not None

    def store_entry(self, record):
        session = self.session
        try:
            for entry in record:
                new_relationship = GOTermRelationship(
                    go_term_1_id=entry['go_term_1_id'],
                    go_term_2_id=entry['go_term_2_id'],
                    information_content_1=entry['information_content_1'],
                    information_content_2=entry['information_content_2'],
                    resnik_distance=entry['resnik_distance'],
                    minimum_branch_length=entry['minimum_branch_length']
                )
                session.add(new_relationship)
            session.commit()
        except Exception as e:
            self.logger.error(f"Error saving GO term relationship: {e}")
            session.rollback()
