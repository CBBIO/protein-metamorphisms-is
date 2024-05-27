import networkx as nx
from sqlalchemy import text
from sqlalchemy.orm import sessionmaker
import obonet

from protein_metamorphisms_is.operations.base.operator import OperatorBase
from protein_metamorphisms_is.sql.model import (
    Protein,
    ProteinGOTermAssociation,
    SequenceGOPrediction,
    GOPerProteinSemanticDistance,
    PredictionMethod,
    EmbeddingType, Sequence
)

class GoPredictionMetricsPerProtein(OperatorBase):
    def __init__(self, conf):
        super().__init__(conf)
        self.logger.info("GoPredictionMetricsPerProtein instance created")
        self.go_graph = obonet.read_obo(self.conf['obo'])
        self.k = conf.get('k', 5)

    def start(self):
        try:
            self.process_protein_metrics()
        except Exception as e:
            self.logger.error(f"Error during GO prediction metrics processing: {e}")
            raise

    def find_descendants(self, nodes):
        descendants = set()
        for node in nodes:
            if node in self.go_graph:
                current_descendants = nx.descendants(self.go_graph, node)
                descendants.update(current_descendants)
        return descendants

    def create_subgraph(self, nodes):
        descendants = self.find_descendants(nodes)
        subgraph_nodes = nodes.union(descendants)
        return self.go_graph.subgraph(subgraph_nodes)

    def process_protein_metrics(self):
        proteins = (self.session.query(Protein)
                    .join(SequenceGOPrediction, Protein.sequence_id == SequenceGOPrediction.sequence_id)
                    .distinct().all())
        prediction_methods = self.session.query(PredictionMethod).all()
        embedding_types = self.session.query(EmbeddingType).all()

        for protein in proteins:
            original = {assoc.go_term.go_id for assoc in protein.go_term_associations}

            for method in prediction_methods:
                for embedding in embedding_types:
                    predictions = {pred.go_term.go_id for pred in self.session.query(SequenceGOPrediction).filter_by(
                        sequence_id=protein.sequence_id,
                        prediction_method_id=method.id,
                        embedding_type_id=embedding.id
                    )}

                    if predictions:
                        original_subgraph_nodes = original.union(self.find_descendants(original))
                        predicted_subgraph_nodes = predictions.union(self.find_descendants(predictions))

                        # Calcula la distancia Jaccard usando los conjuntos de nodos
                        jaccard_distance = self.jaccard_distance(original_subgraph_nodes, predicted_subgraph_nodes)

                        # Guardar la distancia Jaccard
                        self.save_jaccard_distance(protein.entry_name, jaccard_distance, embedding.id, method.id)

    def jaccard_distance(self, set1, set2):
        intersection = len(set1.intersection(set2))
        union = len(set1.union(set2))
        return 1 - intersection / union if union != 0 else 0

    def save_jaccard_distance(self, protein_entry_name, jaccard_distance, embedding_type_id, prediction_method_id):

        distance_record = GOPerProteinSemanticDistance(
            protein_entry_name=protein_entry_name,
            group_distance=jaccard_distance,
            embedding_type_id=embedding_type_id,
            prediction_method_id=prediction_method_id
        )
        self.session.add(distance_record)
        self.session.commit()

