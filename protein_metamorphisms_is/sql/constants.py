# from protein_metamorphisms_is.sql.model.model import StructuralComplexityLevel, StructuralAlignmentType,
# SequenceEmbeddingType, StructureEmbeddingType, PredictionMethod
#
#
#
# def handle_structural_complexity_levels(session, constants):
#     structural_complexity_levels = constants['structural_complexity_levels']
#
#     for level_data in structural_complexity_levels:
#         exists = session.query(StructuralComplexityLevel).filter_by(name=level_data['name']).first()
#         if not exists:
#             complexity_level = StructuralComplexityLevel(**level_data)
#             session.add(complexity_level)
#     session.commit()
#
#
from protein_metamorphisms_is.sql.model.entities.embedding.sequence_embedding import SequenceEmbeddingType
from protein_metamorphisms_is.sql.model.operational.structural_alignment.structural_alignment_type import \
    StructuralAlignmentType


def handle_structural_alignment_types(session, constants):
    structural_alignment_types = constants['structural_alignment_types']

    for level_data in structural_alignment_types:
        exists = session.query(StructuralAlignmentType).filter_by(name=level_data['name']).first()
        if not exists:
            alignment_type = StructuralAlignmentType(**level_data)
            session.add(alignment_type)
    session.commit()

#
#


def handle_sequence_embedding_types(session, constants):
    sequence_embedding_types = constants['sequence_embedding_types']

    for type_data in sequence_embedding_types:
        exists = session.query(SequenceEmbeddingType).filter_by(name=type_data['name']).first()
        if not exists:
            embedding_type = SequenceEmbeddingType(**type_data)
            session.add(embedding_type)
    session.commit()

#
# def handle_structure_embedding_types(session, constants):
#     structure_embedding_types = constants['structure_embedding_types']
#
#     for type_data in structure_embedding_types:
#         exists = session.query(StructureEmbeddingType).filter_by(name=type_data['name']).first()
#         if not exists:
#             embedding_type = StructureEmbeddingType(**type_data)
#             session.add(embedding_type)
#     session.commit()
#
#
# def handle_prediction_methods(session, constants):
#     prediction_methods = constants['prediction_methods']
#
#     for method_data in prediction_methods:
#         exists = session.query(PredictionMethod).filter_by(name=method_data['name']).first()
#         if not exists:
#             prediction_method = PredictionMethod(**method_data)
#             session.add(prediction_method)
#     session.commit()
