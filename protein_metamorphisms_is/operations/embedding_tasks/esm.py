from transformers import AutoTokenizer, EsmModel
import torch


def embedding_task(session,chains,module,model_name):
    # Verificar si CUDA está disponible
    if not torch.cuda.is_available():
        raise Exception("CUDA is not available. This script requires a GPU with CUDA.")

    # Configurar el dispositivo
    device = torch.device("cuda")

    # Cargar el tokenizador y el modelo
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = EsmModel.from_pretrained(model_name).to(device)

    # Preparar la secuencia

    with torch.no_grad():  # Desactivar el cálculo de gradientes
        for chain in chains:
            tokens = tokenizer(chain.sequence, return_tensors="pt", truncation=True, padding=True)
        # Mover los tokens al dispositivo correcto
            tokens = {k: v.to(device) for k, v in tokens.items()}

            # Obtener los embeddings del modelo
            outputs = model(**tokens)
            embeddings = outputs.last_hidden_state
            # embeddings es un tensor de shape (batch_size, sequence_length, hidden_size)

            print(embeddings.shape)

