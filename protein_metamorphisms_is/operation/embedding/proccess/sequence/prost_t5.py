import traceback
from retry import retry
from transformers import T5Tokenizer, T5EncoderModel
import re
import torch

def load_model(model_name):
    return T5EncoderModel.from_pretrained(model_name).to(torch.device("cuda"))

def load_tokenizer(model_name):
    return T5Tokenizer.from_pretrained(model_name, do_lower_case=False)

def embedding_task(sequences, model, tokenizer, batch_size=32):
    if not torch.cuda.is_available():
        raise Exception("CUDA is not available. This script requires a GPU with CUDA.")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.eval()
    embedding_records = []

    # Preprocess sequences
    sequences_processed = [
        ("<AA2fold> " + " ".join(list(re.sub(r"[UZOB]", "X", seq))) if seq.isupper()
         else "<fold2AA> " + " ".join(list(seq)))
        for seq in sequences
    ]

    # Process sequences in batches
    for i in range(0, len(sequences_processed), batch_size):
        batch_sequences = sequences_processed[i:i + batch_size]
        inputs = tokenizer.batch_encode_plus(
            batch_sequences,
            padding="longest",
            truncation=True,
            max_length=512,
            add_special_tokens=True,
            return_tensors="pt"
        ).to(device)

        with torch.no_grad():
            try:
                outputs = model(input_ids=inputs.input_ids, attention_mask=inputs.attention_mask)
                embeddings = outputs.last_hidden_state.mean(dim=1)

                # Collect embeddings for the batch
                for idx, sequence in enumerate(batch_sequences):
                    record = {
                        'sequence': sequences[i + idx],  # Original sequence
                        'embedding': embeddings[idx].cpu().numpy().tolist(),  # Embedding vector
                        'shape': embeddings[idx].shape
                    }
                    embedding_records.append(record)

            except Exception as e:
                print(f"Error processing batch {i // batch_size}: {e}")
                continue

    return embedding_records
