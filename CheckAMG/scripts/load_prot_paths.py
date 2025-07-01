#!/usr/bin/env python3

import os

def load_prots(protein_dir):
    vmag_proteins_subdir = os.path.join(protein_dir, 'vMAG_proteins')
    
    if os.path.exists(os.path.join(protein_dir, 'single_contig_proteins.faa')):
        if os.path.exists(vmag_proteins_subdir) and os.path.isdir(vmag_proteins_subdir):
            prots = [os.path.join(protein_dir, 'single_contig_proteins.faa')] + [os.path.join(vmag_proteins_subdir, f) for f in os.listdir(vmag_proteins_subdir) if f.endswith('.faa')]
        else:
            prots = [os.path.join(protein_dir, 'single_contig_proteins.faa')]
    else:
        if os.path.exists(vmag_proteins_subdir) and os.path.isdir(vmag_proteins_subdir):
            prots = [os.path.join(vmag_proteins_subdir, f) for f in os.listdir(vmag_proteins_subdir) if f.endswith('.faa')]
        else:
            raise FileNotFoundError(f"No valid protein files found in search path {protein_dir} or {vmag_proteins_subdir}.")
        
    return prots