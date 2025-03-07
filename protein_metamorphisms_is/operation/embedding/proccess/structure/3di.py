from mini3di import Encoder


def load_model(model_name, conf):
    # mini3di doesn't require loading a model in the traditional sense
    return Encoder()


def load_tokenizer(model_name):
    # mini3di doesn't use a tokenizer, but we provide a dummy function to maintain interface consistency
    return None


def embedding_task(structures, model, tokenizer):
    """
    Processes structures to generate embeddings using mini3di.

    Args:
    - structures: List of structures to process.
    - model: Preloaded mini3di Encoder.
    - tokenizer: Not used in mini3di, kept for interface compatibility.

    Returns:
    A list of embedding records.
    """
    embedding_records = []
    for structure in structures:
        # Assume structure is a dictionary with the required atomic coordinates
        try:
            ca = structure['ca']
            cb = structure['cb']
            n = structure['n']
            c = structure['c']

            states = model.encode_atoms(ca=ca, cb=cb, n=n, c=c)
            embedding_shape = states.shape

            # Prepare the record
            record = {
                'structure': structure,
                'embedding': states.tolist(),
                'shape': embedding_shape
            }

            embedding_records.append(record)
        except Exception as e:
            print(f"Failed to process structure {structure}: {e}")
            continue

    return embedding_records
