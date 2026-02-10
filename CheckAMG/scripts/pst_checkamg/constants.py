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
