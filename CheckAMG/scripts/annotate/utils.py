# CheckAMG/scripts/common/utils.py

from __future__ import annotations
import os
from pathlib import Path


def load_prots(protein_dir, bin_proteins_subdir):    
    if os.path.exists(os.path.join(protein_dir, 'single_contig_proteins.faa')):
        if os.path.exists(bin_proteins_subdir) and os.path.isdir(bin_proteins_subdir):
            prots = [os.path.join(protein_dir, 'single_contig_proteins.faa')] + [os.path.join(bin_proteins_subdir, f) for f in os.listdir(bin_proteins_subdir) if f.endswith('.faa')]
        else:
            prots = [os.path.join(protein_dir, 'single_contig_proteins.faa')]
    else:
        if os.path.exists(bin_proteins_subdir) and os.path.isdir(bin_proteins_subdir):
            prots = [os.path.join(bin_proteins_subdir, f) for f in os.listdir(bin_proteins_subdir) if f.endswith('.faa')]
        else:
            raise FileNotFoundError(f"No valid protein files found in search path {protein_dir} or {bin_proteins_subdir}.")
        
    return prots


def as_parquet_path_if_enabled(path: str | Path, save_to_parquet: bool) -> Path:
    p = Path(path)
    if not save_to_parquet:
        return p

    # already parquet
    if p.suffix.lower() == ".parquet":
        return p

    # common expected case: tsv -> parquet
    if p.suffix.lower() == ".tsv":
        return p.with_suffix(".parquet")

    # if it's some other suffix, replace it
    if p.suffix:
        return p.with_suffix(".parquet")

    # no suffix at all
    return p.with_suffix(".parquet")