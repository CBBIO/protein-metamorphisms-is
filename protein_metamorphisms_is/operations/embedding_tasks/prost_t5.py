from transformers import T5Tokenizer, T5EncoderModel
import re
import torch

from protein_metamorphisms_is.sql.model import ChainEmbedding


def embedding_task(session,chains,model_name, embedding_type_id):
    if not torch.cuda.is_available():
        raise Exception("CUDA is not available. This script requires a GPU with CUDA.")

    device = torch.device("cuda")
    tokenizer = T5Tokenizer.from_pretrained(model_name, do_lower_case=False)
    model = T5EncoderModel.from_pretrained(model_name).to(device)
    model.eval()

    with torch.no_grad():
        for chain in chains:
            sequence_processed = " ".join(list(re.sub(r"[UZOB]", "X", chain.sequence)))
            sequence_processed = "<AA2fold> " + sequence_processed if sequence_processed.isupper() else "<fold2AA> " + sequence_processed
            inputs = tokenizer(sequence_processed, return_tensors="pt", padding=True, truncation=True, max_length=512, add_special_tokens=True).to(device)

            outputs = model(input_ids=inputs.input_ids, attention_mask=inputs.attention_mask)
            embeddings = outputs.last_hidden_state.mean(dim=1).cpu().numpy().tolist()[0]

            embedding_entry = ChainEmbedding(
                pdb_chain_id=chain.id,
                embedding_type_id=embedding_type_id,
                embedding=embeddings
            )
            session.add(embedding_entry)

    session.commit()