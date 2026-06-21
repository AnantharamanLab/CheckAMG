import os
os.environ["NUMEXPR_MAX_THREADS"] = str(snakemake.threads)

import tables as tb

from CheckAMG.scripts.denovo.data import GenomeDataModule
from CheckAMG.scripts.denovo.embed import embed
from CheckAMG.scripts.denovo.utils import (
    find_best_checkpoint,
    load_model,
    _extract_label_from_datamodule,
)

from CheckAMG.scripts.common.snakemake_logging import init_snakemake_logger, emit_step_banner
from CheckAMG.scripts.common.mem_limit import set_memory_limit

log_file = snakemake.params.log
logger = init_snakemake_logger(snakemake, redirect_streams=False)
emit_step_banner(title=f"Step {snakemake.params.step[0]}/{snakemake.params.step[1]}: Embed training proteins with the trained CheckAMG-PST model", log_path=log_file, append_raw=False)
set_memory_limit(snakemake.resources.mem)
logger.debug(f"Memory limit set to {snakemake.resources.mem} GB")

import warnings
from tables import UnclosedFileWarning
warnings.filterwarnings("ignore", category=UnclosedFileWarning, module="tables")


def main():
    train_file = snakemake.params.train_file
    train_embed_output = snakemake.params.train_embed_output
    logdir = snakemake.params.logdir

    batch_size = snakemake.params.batch_size
    class_weighting = snakemake.params.class_weighting
    accelerator = snakemake.params.accelerator

    ckpt = find_best_checkpoint(logdir)
    logger.info(f"Loading best trained model from checkpoint: {ckpt}")
    model = load_model(ckpt, device=accelerator)

    datamodule = GenomeDataModule(train_file=train_file, class_weighting=class_weighting)
    datamodule.setup()

    logger.info(f"Embedding training proteins -> {train_embed_output}")
    embed(
        model=model,
        dataset=datamodule.train_dataset,  # type: ignore
        output=train_embed_output,
        batch_size=batch_size,
    )

    # Check if the labels are already present in the output PST embedding file, and if not, add them
    label_exists = False
    with tb.open_file(train_embed_output, mode="r") as fp:
        label_exists = "label" in fp.root

    if not label_exists:
        # Add the labels from the ESM embeddings to the PST embeddings, so they can be used for inference
        logger.debug(f"Adding labels from the ESM embeddings to the PST embeddings...")
        train_labels = _extract_label_from_datamodule(datamodule)
        with tb.open_file(train_embed_output, mode="a") as fp:
            atom = tb.Atom.from_dtype(train_labels.numpy().dtype)
            ds = fp.create_carray(fp.root, "label", atom, train_labels.shape)
            ds[:] = train_labels.numpy()
    else:
        logger.debug(f"Labels already exist in the PST embedding file. Skipping label addition.")

    datamodule.close("train")


if __name__ == "__main__":
    main()
