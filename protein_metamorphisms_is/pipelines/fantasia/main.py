
from protein_metamorphisms_is.pipelines.fantasia.embedder import SequenceEmbedder
from protein_metamorphisms_is.pipelines.fantasia.lookup import EmbeddingLookUp
from protein_metamorphisms_is.helpers.config.yaml import read_yaml_config
from protein_metamorphisms_is.operation.extraction.accessions import AccessionManager
from protein_metamorphisms_is.operation.extraction.uniprot import UniProtExtractor
from protein_metamorphisms_is.operation.embedding.sequence_embedding import SequenceEmbeddingManager

from protein_metamorphisms_is.helpers.config.yaml import read_yaml_config

from protein_metamorphisms_is.operation.extraction.pdb import PDBExtractor
from protein_metamorphisms_is.operation.embedding.structure_3di import Structure3DiManager
from protein_metamorphisms_is.operation.clustering.sequence_clustering import SequenceClustering
from protein_metamorphisms_is.operation.clustering.structural_subclustering import StructuralSubClustering
from protein_metamorphisms_is.operation.functional.annotation_transfer.sequence_go_annotation import \
    SequenceGOAnnotation
from protein_metamorphisms_is.operation.structural_alignment.structural_alignment import StructuralAlignmentManager


from protein_metamorphisms_is.operation.functional.multifunctionality.go_multifunctionality_metrics import \
    GoMultifunctionalityMetrics



def main(config_path="/home/bioxaxi/PycharmProjects/protein-metamorphisms-is/protein_metamorphisms_is/pipelines/fantasia/config.yaml"):
    conf = read_yaml_config(config_path)
    # AccessionManager(conf).load_accessions_from_csv()
    # UniProtExtractor(conf).start()
    # SequenceEmbeddingManager(conf).start()
    embedder = SequenceEmbedder(conf)
    embedder.start()
    lookup = EmbeddingLookUp(conf)
    lookup.start()





if __name__ == "__main__":
    main()
