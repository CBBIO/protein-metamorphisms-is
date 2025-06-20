from transformers import T5Tokenizer, T5EncoderModel
import torch
import numpy as np
import gc


def load_model(model_name, conf):
    device = torch.device(conf['embedding'].get('device', "cuda"))
    dtype = torch.float16 if device.type == "cuda" else torch.float32
    model = T5EncoderModel.from_pretrained(model_name, torch_dtype=dtype).to(device).eval()
    return model


def load_tokenizer(model_name):
    return T5Tokenizer.from_pretrained(model_name)


def embedding_task(sequences, model, tokenizer, device, batch_size=8, embedding_type_id=None):
    model.to(device)
    model.eval()
    embedding_records = []

    for i in range(0, len(sequences), batch_size):
        batch = sequences[i:i + batch_size]
        processed = ["[NLU]" + seq["sequence"] for seq in batch]

        inputs = tokenizer.batch_encode_plus(
            processed,
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=512
        ).to(device)

        with torch.no_grad():
            try:
                outputs = model(**inputs)
                embeddings = outputs.last_hidden_state.mean(dim=1).cpu().numpy()

                for idx, seq in enumerate(batch):
                    embedding_records.append({
                        "sequence_id": seq["sequence_id"],
                        "embedding_type_id": embedding_type_id,
                        "sequence": seq["sequence"],
                        "embedding": embeddings[idx].astype(np.float16),
                        "shape": embeddings[idx].shape
                    })
            except Exception as e:
                print(f"⚠️ Error en batch {i // batch_size}: {e}")
                continue

        torch.cuda.empty_cache()
        gc.collect()

    return embedding_records
