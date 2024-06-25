from protein_metamorphisms_is.helpers.config.yaml import read_yaml_config
from protein_metamorphisms_is.information_system.pdb import PDBExtractor
from protein_metamorphisms_is.information_system.uniprot import UniProtExtractor
from protein_metamorphisms_is.operations.cdhit import CDHit
from protein_metamorphisms_is.operations.go_prediction import GoPrediction
from protein_metamorphisms_is.operations.protein_go_prediction_metrics import GoPredictionMetricsPerProtein

from protein_metamorphisms_is.operations.seq_embeddings import SequenceEmbeddingManager
from protein_metamorphisms_is.operations.optics import OpticsClustering

from protein_metamorphisms_is.operations.go_metrics import GoMetrics
from protein_metamorphisms_is.operations.seq_embeddings import SequenceEmbeddingManager

from protein_metamorphisms_is.operations.structural_alignment import StructuralAlignmentManager


def main(config_path="config/config.yaml"):
    conf = read_yaml_config(config_path)
    UniProtExtractor(conf).start()
    GoMetrics(conf).start()
    # PDBExtractor(conf).start()
    # SequenceEmbeddingManager(conf).start()
    # CDHit(conf).start()
    # OpticsClustering(conf).start()
    # StructuralAlignmentManager(conf).start()
    # GoPrediction(conf).start()
    # GoPredictionMetricsPerProtein(conf).start()


if __name__ == "__main__":
    main()
