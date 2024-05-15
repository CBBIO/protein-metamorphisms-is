from transformers import AutoTokenizer, EsmModel
import torch

from protein_metamorphisms_is.sql.model import SequenceEmbedding


def embedding_task(session, sequences, model_name, embedding_type_id):
    if not torch.cuda.is_available():
        raise Exception("CUDA is not available. This script requires a GPU with CUDA.")

    device = torch.device("cuda")
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = EsmModel.from_pretrained(model_name).to(device)

    with torch.no_grad():
        for sequence in sequences:
            tokens = tokenizer(sequence.sequence, return_tensors="pt", truncation=True, padding=True)
            tokens = {k: v.to(device) for k, v in tokens.items()}
            outputs = model(**tokens)
            embeddings = outputs.last_hidden_state.mean(dim=1)
            embedding_shape = embeddings.shape

            embedding_entry = SequenceEmbedding(
                sequence_id=sequence.id,
                embedding_type_id=embedding_type_id,
                embedding=embeddings.cpu().numpy().tolist()[0],
                shape=embedding_shape
            )
            session.add(embedding_entry)

    session.commit()
