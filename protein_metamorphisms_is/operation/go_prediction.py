from sqlalchemy import text

from sqlalchemy.orm import Session

from protein_metamorphisms_is.operation.base.operator import OperatorBase
from protein_metamorphisms_is.sql.model import SequenceEmbedding, ProteinGOTermAssociation, SequenceGOPrediction, PredictionMethod


class GoPrediction(OperatorBase):
    def __init__(self, conf):
        super().__init__(conf)
        self.logger.info("GoPrediction instance created")
        self.k = conf.get('k', 5)  # Permitir configuración externa de k

    def start(self):
        try:
            embeddings = self.session.query(SequenceEmbedding).all()
            prediction_method = self.session.query(PredictionMethod).filter_by(name="Cosine Similarity").first()

            for embedding in embeddings:
                self.logger.info(f"Processing embedding ID {embedding.id}")
                nearest_sequence_ids = self.find_k_nearest_sequences(embedding.id, embedding.embedding_type_id)
                predicted_go_terms_info = self.fetch_associated_go_terms(nearest_sequence_ids)

                if predicted_go_terms_info:
                    self.save_go_predictions(embedding.sequence_id, predicted_go_terms_info, embedding.embedding_type_id, prediction_method.id)
                    self.logger.info(f"Saved GO terms for embedding ID {embedding.id}")
                else:
                    self.logger.info(f"No GO terms found for embedding ID {embedding.id}")
        except Exception as e:
            self.logger.error(f"Error during GO prediction process: {e}")
            raise

    def find_k_nearest_sequences(self, sequence_embedding_id, embedding_type_id):
        query = text("""
        SELECT sequence_id
        FROM sequence_embeddings
        WHERE id != :seq_id AND embedding_type_id = :emb_type_id
        ORDER BY embedding <-> (
            SELECT embedding
            FROM sequence_embeddings
            WHERE id = :seq_id AND embedding_type_id = :emb_type_id
        ) LIMIT :k;
        """)
        result = self.session.execute(query, {'seq_id': sequence_embedding_id, 'emb_type_id': embedding_type_id, 'k': self.k}).fetchall()
        return [r[0] for r in result]

    def fetch_associated_go_terms(self, sequence_ids):
        if not sequence_ids:
            return []

        query = text("""
        SELECT
            s.id AS sequence_id,
            p.entry_name AS protein_entry_name,
            gt.go_id
        FROM
            proteins p
        JOIN
            sequences s ON p.sequence_id = s.id
        JOIN
            protein_go_term_association pgta ON p.entry_name = pgta.protein_entry_name
        JOIN
            go_terms gt ON pgta.go_id = gt.go_id
        WHERE
            s.id IN :seq_ids;
        """)
        result = self.session.execute(query, {'seq_ids': tuple(sequence_ids)}).fetchall()
        return result

    def save_go_predictions(self, sequence_id, predicted_go_terms_info, embedding_type_id, prediction_method_id):
        for ref_sequence_id, ref_protein_entry_name, go_id in predicted_go_terms_info:
            prediction = self.session.query(SequenceGOPrediction).filter_by(
                sequence_id=sequence_id,
                ref_protein_entry_name=ref_protein_entry_name,
                go_id=go_id,
                k=self.k,
                embedding_type_id=embedding_type_id,
                prediction_method_id=prediction_method_id
            ).one_or_none()

            if prediction:
                self.logger.info(f"Update existing entry for {ref_protein_entry_name}")
                # Actualizar campos si es necesario
            else:
                self.logger.info(f"Adding new prediction for {ref_protein_entry_name}")
                prediction = SequenceGOPrediction(
                    sequence_id=sequence_id,
                    ref_protein_entry_name=ref_protein_entry_name,
                    go_id=go_id,
                    k=self.k,
                    embedding_type_id=embedding_type_id,
                    prediction_method_id=prediction_method_id
                )
                self.session.add(prediction)

            try:
                self.session.commit()
            except Exception as e:
                self.logger.error(f"Error during commit: {e}")
                self.session.rollback()
