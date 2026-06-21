from pst.typing import EdgeIndexStrategy

NUM_CLASSES = 4
TRUE_AVG_IDX = 3
EMBED_FIELD = "ctx_ptn"
LABEL_FIELD = "label"

# PST-TL-P__large datamodule args
CHUNK_SIZE = 25
FRAGMENT_SIZE = -1
THRESHOLD = -1
EDGE_STRATEGY = EdgeIndexStrategy.chunked

# ESM2 embedding resilience: the filtered FASTA is embedded in scaffold-aware
# chunks of this many proteins so that a transient GPU/CUDA fault only loses the
# in-progress chunk (completed chunks are skipped on re-run). Each chunk is
# embedded then graphified independently and merged into the final graph file.
ESM2_EMBED_CHUNK_SIZE = 250_000
# Number of attempts per chunk before giving up (in-process retry on CUDA errors).
ESM2_EMBED_RETRIES = 3
