from transformers import T5Tokenizer, T5EncoderModel
import torch
import re

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

# Load the tokenizer
tokenizer = T5Tokenizer.from_pretrained('Rostlab/ProstT5', do_lower_case=False)

# Load the model
model = T5EncoderModel.from_pretrained("Rostlab/ProstT5").to(device)

# only GPUs support half-precision currently; if you want to run on CPU use full-precision (not recommended, much slower)
if device.type == 'cpu':
    model.float()
else:
    model.half()

# Prepare your protein sequences/structures as a list.
# Replace all rare/ambiguous amino acids by X (3Di sequences do not have those) and introduce white-space between all sequences (AAs and 3Di)
sequence_examples = ["PRTEINO", "strct", "SEQUENCE1", "sequence2", "SEQ3", "seq4", "SEQUENCE5", "seq6"]
sequence_examples = [" ".join(list(re.sub(r"[UZOB]", "X", sequence))) for sequence in sequence_examples]

# The direction of the translation is indicated by two special tokens:
# if you go from AAs to 3Di (or if you want to embed AAs), you need to prepend "<AA2fold>"
# if you go from 3Di to AAs (or if you want to embed 3Di), you need to prepend "<fold2AA>"
sequence_examples = ["<AA2fold> " + s if s.isupper() else "<fold2AA> " + s for s in sequence_examples]

# Function to process each sequence individually
def process_sequence(sequence):
    # Tokenize sequence and pad up to the longest sequence in the batch
    ids = tokenizer.encode_plus(sequence,
                                add_special_tokens=True,
                                padding="longest",
                                return_tensors='pt').to(device)

    # Generate embeddings
    with torch.no_grad():
        embedding_repr = model(
            ids.input_ids,
            attention_mask=ids.attention_mask
        )

    # Extract residue embeddings for the sequence and remove padded & special tokens, incl. prefix
    emb = embedding_repr.last_hidden_state[0, 1:1 + len(sequence.split()) - 1]  # adjust slicing as needed
    return emb

# Process each sequence individually and store results
individual_embeddings = []
for seq in sequence_examples:
    embedding = process_sequence(seq)
    individual_embeddings.append(embedding)
    print(f"Individual embedding for sequence '{seq}':")
    print(embedding)

# Tokenize sequences in batch and pad up to the longest sequence in the batch
batch_ids = tokenizer.batch_encode_plus(sequence_examples,
                                        add_special_tokens=True,
                                        padding="longest",
                                        return_tensors='pt').to(device)

# Generate batch embeddings
with torch.no_grad():
    batch_embedding_repr = model(
        batch_ids.input_ids,
        attention_mask=batch_ids.attention_mask
    )

# Extract residue embeddings for each sequence in the batch and remove padded & special tokens, incl. prefix
batch_embeddings = []
for i, seq in enumerate(sequence_examples):
    emb = batch_embedding_repr.last_hidden_state[i, 1:1 + len(seq.split()) - 1]  # adjust slicing as needed
    batch_embeddings.append(emb)
    print(f"Batch embedding for sequence '{seq}':")
    print(emb)

# Compare individual and batch embeddings
for i, seq in enumerate(sequence_examples):
    # Calculate mean embeddings per protein
    emb_individual_mean = individual_embeddings[i].mean(dim=0)  # shape (1024)
    emb_batch_mean = batch_embeddings[i].mean(dim=0)  # shape (1024)

    # Calculate euclidean distance
    distance = torch.dist(emb_individual_mean, emb_batch_mean)

    print(f"Comparison for sequence '{seq}':")
    print(f"Mean individual embedding: {emb_individual_mean}")
    print(f"Mean batch embedding: {emb_batch_mean}")
    print(f"Euclidean distance between embeddings: {distance.item()}")

# Calculate distance between sequence pairs
for i in range(len(sequence_examples)):
    for j in range(i+1, len(sequence_examples)):
        indiv_dist = torch.dist(individual_embeddings[i].mean(dim=0), individual_embeddings[j].mean(dim=0))
        batch_dist = torch.dist(batch_embeddings[i].mean(dim=0), batch_embeddings[j].mean(dim=0))
        print(f"Distance between sequence {i} and {j} in individual embeddings: {indiv_dist.item()}")
        print(f"Distance between sequence {i} and {j} in batch embeddings: {batch_dist.item()}")
