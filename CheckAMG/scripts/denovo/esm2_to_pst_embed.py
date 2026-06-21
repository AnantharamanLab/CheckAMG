from pathlib import Path

import os
os.environ["NUMEXPR_MAX_THREADS"] = str(snakemake.threads)

import faiss
import pst
import torch

from CheckAMG.scripts.denovo.constants import (
    CHUNK_SIZE,
    EDGE_STRATEGY,
    FRAGMENT_SIZE,
    THRESHOLD,
)
from CheckAMG.scripts.denovo.embed import embed
from CheckAMG.scripts.denovo.model import CheckAMGPST
from CheckAMG.scripts.denovo.utils import load_model

from CheckAMG.scripts.common.snakemake_logging import init_snakemake_logger, emit_step_banner
from CheckAMG.scripts.common.mem_limit import set_memory_limit

log_file = snakemake.params.log
logger = init_snakemake_logger(snakemake, redirect_streams=False)
emit_step_banner(title=f"Step {snakemake.params.step[0]}/{snakemake.params.step[1]}: Compute PST embeddings from ESM2", log_path=log_file, append_raw=False)
set_memory_limit(snakemake.resources.mem)
logger.debug(f"Memory limit set to {snakemake.resources.mem} GB")

import warnings
from tables import UnclosedFileWarning
warnings.filterwarnings("ignore", category=UnclosedFileWarning, module="tables")

OptPath = None | Path

def compute_pst_embeddings_from_esm2(
    model: CheckAMGPST,
    query_file_esm2: OptPath,
    query_file_pst: OptPath,
    batch_size: int,
):
    
    logger.info(f"Building a PST LazyGenomeDataset from ESM2 graph file {query_file_esm2}...")
    dataset = pst.LazyGenomeDataset(
            file=query_file_esm2,
            edge_strategy=EDGE_STRATEGY,
            chunk_size=CHUNK_SIZE,
            threshold=THRESHOLD,
            fragment_size=FRAGMENT_SIZE,
            )

    try:
        logger.info(f"Computing PST embeddings for {query_file_esm2}...")
        embed(
            model=model,
            dataset=dataset,
            output=query_file_pst,
            batch_size=batch_size,
        )
    except Exception as e:
        logger.error(f"Failed to compute PST embeddings for {query_file_esm2}: {e}")
        raise RuntimeError(f"Failed to compute PST embeddings for {query_file_esm2}") from e

    logger.info(f"PST embeddings computed and saved to {query_file_pst}.")

def _handle_cpu_accelerator(accelerator: str | None, devices: int):
    if accelerator == "cpu" or (
        accelerator == "auto" and not torch.cuda.is_available()
    ):
        torch.set_num_threads(devices)
        faiss.omp_set_num_threads(devices)

def main():
    logger.info("PST embedding starting...")

    accelerator = str(snakemake.params.accelerator)
    devices = int(snakemake.params.devices)
    logger.debug(f"Accelerator={accelerator}, devices={devices}")
    _handle_cpu_accelerator(accelerator, devices)
    
    model_ckpt = snakemake.params.model_ckpt
    logger.debug(f"Model checkpoint for PST embedding: {model_ckpt}")
    model = load_model(ckpt=model_ckpt, device=accelerator)
    
    esm2_graph_file = Path(snakemake.params.esm2_graph_file)
    query_file_pst = Path(snakemake.params.query_file_pst)
    logger.debug(f"ESM2 graph file for PST embedding: {esm2_graph_file}, Output PST embedding file: {query_file_pst}")

    compute_pst_embeddings_from_esm2(
        model=model,
        query_file_esm2=esm2_graph_file,
        query_file_pst=query_file_pst,
        batch_size=snakemake.params.batch_size
        )
    
    logger.info("PST embedding completed.")

if __name__ == "__main__":
    main()