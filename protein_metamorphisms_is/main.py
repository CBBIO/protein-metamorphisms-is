from protein_metamorphisms_is.helpers.config.yaml import read_yaml_config

from protein_metamorphisms_is.operation.extraction.accessions import AccessionManager
from protein_metamorphisms_is.operation.extraction.pdb import PDBExtractor
from protein_metamorphisms_is.operation.extraction.uniprot import UniProtExtractor
from protein_metamorphisms_is.operation.embedding.sequence_embedding import SequenceEmbeddingManager
from protein_metamorphisms_is.operation.embedding.structure_3di import Structure3DiManager
from protein_metamorphisms_is.operation.clustering.sequence_clustering import SequenceClustering
from protein_metamorphisms_is.operation.clustering.structural_subclustering import StructuralSubClustering
from protein_metamorphisms_is.operation.functional.annotation_transfer.sequence_go_annotation import \
    SequenceGOAnnotation
from protein_metamorphisms_is.operation.structural_alignment.structural_alignment import StructuralAlignmentManager


from protein_metamorphisms_is.operation.functional.multifunctionality.go_multifunctionality_metrics import \
    GoMultifunctionalityMetrics

def main(config_path="config/config.yaml"):
    conf = read_yaml_config(config_path)
    AccessionManager(conf).fetch_accessions_from_api()
    AccessionManager(conf).load_accessions_from_csv()
    UniProtExtractor(conf).start()
    PDBExtractor(conf).start()
    SequenceEmbeddingManager(conf).start()
    Structure3DiManager(conf).start()
    SequenceClustering(conf).start()
    StructuralSubClustering(conf).start()
    StructuralAlignmentManager(conf).start()
    SequenceGOAnnotation(conf).start()
    # GoMultifunctionalityMetrics(conf).start()




if __name__ == "__main__":
    main()
