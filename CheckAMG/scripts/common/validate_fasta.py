# CheckAMG/scripts/common/validate_fasta.py

import re
from CheckAMG.scripts.common.fasta_io import iter_fasta_header_lines

# Prodigal protein headers separate fields with whitespace-padded ' # '.
_PRODIGAL_HEADER_SPLIT_RE = re.compile(r"\s+#\s+")


def is_prodigal_header(header_text):
    """
    True iff `header_text` (FASTA header content after '>') matches Prodigal's
    protein header format:
        seq_id # start # end # strand # ID=...
    The seq_id must be a single whitespace-free token whose final
    underscore-delimited segment is a non-negative integer (the ORF index).
    Examples accepted:
        viral_genome_1115664__k141_540481_123 # 1 # 500 # -1 # ID=...
        viral_genome_1115664_k141_540481_1_123 # 1 # 500 # -1 # ID=...
    Examples rejected:
        viral_genome_1115664 k141_540481_123 # 1 # 500 # -1 # ID=...
        viral_genome_1115664_k141_540481 # 1 # 500 # -1 # ID=...   (no ORF index)
    """
    parts = _PRODIGAL_HEADER_SPLIT_RE.split(header_text.strip())
    if len(parts) < 5:
        return False

    seq_id = parts[0]
    if not seq_id or any(ch.isspace() for ch in seq_id):
        return False
    if "_" not in seq_id:
        return False
    if not seq_id.rsplit("_", 1)[1].isdigit():
        return False

    try:
        int(parts[1])
        int(parts[2])
        int(parts[3])
    except ValueError:
        return False

    if "ID=" not in parts[4]:
        return False
    return True


def _prodigal_seq_id(header_text):
    return _PRODIGAL_HEADER_SPLIT_RE.split(header_text.strip(), maxsplit=1)[0]


def assert_unique_fasta_names(path, logger):
    """
    Verify all FASTA records in `path` have unique sequence names. The "name"
    is the first whitespace-delimited token after `>` (the standard FASTA
    convention; also the seq_id portion of a Prodigal header).
    Single header-only pass; sequences are not parsed.
    """
    seen = set()
    n = 0
    for header_text in iter_fasta_header_lines(path):
        stripped = header_text.strip()
        if not stripped:
            continue
        name = stripped.split(None, 1)[0]
        if name in seen:
            logger.error(f"Duplicate FASTA header name detected: {name} in {path}")
            raise RuntimeError(f"Duplicate FASTA header name detected: {name} in {path}")
        seen.add(name)
        n += 1
    if n == 0:
        logger.error(f"FASTA appears empty or has no headers: {path}")
        raise RuntimeError(f"FASTA appears empty or has no headers: {path}")


def validate_fasta_inputs(input_files, input_type, logger, *, label="input"):
    """
    Validate user-supplied FASTA inputs in a single header-only pass per file.

      * 'nucl': each header must be a single token (no internal whitespace).
                Descriptions after the name are not allowed.
      * 'prot': each header must be Prodigal-formatted; the seq_id must be a
                whitespace-free token whose final underscore-delimited
                segment is the ORF index (an integer).

    Also enforces unique sequence names within each file.
    """
    if input_type not in ("nucl", "prot"):
        logger.error(f"Invalid input type: {input_type}. Expected 'nucl' or 'prot'.")
        raise ValueError(f"Invalid input type: {input_type}. Expected 'nucl' or 'prot'.")
    
    seen = {}
    for path in input_files:
        logger.debug(f"Validating FASTA ({input_type}): {path}")

        n = 0
        for header_text in iter_fasta_header_lines(path):
            stripped = header_text.strip()
            if not stripped:
                logger.error(f"{label} FASTA contains an empty header line in {path}")
                raise RuntimeError(f"{label} FASTA contains an empty header line in {path}")

            if input_type == "nucl":
                # if any(ch.isspace() for ch in stripped):
                #     logger.error(
                #         f"{label} FASTA header contains whitespace in sequence name "
                #         f"in {path}: '>{header_text}'"
                #     )
                #     raise RuntimeError(
                #         f"{label} FASTA headers must not contain whitespace in sequence names: {path}"
                #     )
                name = stripped
            else:
                if not is_prodigal_header(stripped):
                    logger.error(
                        f"{label} FASTA header is not in Prodigal format "
                        f"in {path}: '>{header_text}'"
                    )
                    raise RuntimeError(
                        f"{label} FASTA headers must be in Prodigal format "
                        f"(>seq_id # start # end # strand # ID=...): {path}"
                    )
                name = _prodigal_seq_id(stripped)

            if name in seen.keys():
                logger.error(f"Duplicate FASTA header name detected: {name} in {path} and {seen[name]}")
                raise RuntimeError(f"Duplicate FASTA header name detected: {name} in {path} and {seen[name]}")
            seen[name] = str(path)
            n += 1

        if n == 0:
            logger.error(f"{label} FASTA appears empty or has no headers: {path}")
            raise RuntimeError(f"{label} FASTA appears empty or has no headers: {path}")
