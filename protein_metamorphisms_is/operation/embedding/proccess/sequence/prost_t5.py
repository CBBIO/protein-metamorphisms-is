import traceback
from retry import retry
from transformers import T5Tokenizer, T5EncoderModel
import re
import torch

def load_model(model_name):
    return T5EncoderModel.from_pretrained(model_name).to(torch.device("cuda"))

def load_tokenizer(model_name):
    return T5Tokenizer.from_pretrained(model_name, do_lower_case=False)

def embedding_task(sequences, model, tokenizer):
    if not torch.cuda.is_available():
        raise Exception("CUDA is not available. This script requires a GPU with CUDA.")

    device = torch.device("cuda")
    model.eval()
    embedding_records = []

    with torch.no_grad():
        for sequence in sequences:
            sequence_processed = " ".join(list(re.sub(r"[UZOB]", "X", sequence)))
            sequence_processed = "<AA2fold> " + sequence_processed if sequence_processed.isupper() else "<fold2AA> " + sequence_processed
            inputs = tokenizer(sequence_processed, return_tensors="pt", padding=True, truncation=True, max_length=512, add_special_tokens=True).to(device)

            try:
                outputs = model(input_ids=inputs.input_ids, attention_mask=inputs.attention_mask)
                embeddings = outputs.last_hidden_state.mean(dim=1)
                embedding_shape = embeddings.shape

                # Prepare the record
                record = {
                    'sequence': sequence,
                    'embedding': embeddings.cpu().numpy().tolist()[0],
                    'shape': embedding_shape
                }

                embedding_records.append(record)
            except Exception as e:
                traceback.print_exc()  # This will print the traceback
                continue

    return embedding_records
