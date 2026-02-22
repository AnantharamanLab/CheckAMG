#!/usr/bin/env python3

import os

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