from pathlib import Path
from sys import stdout

import os
os.environ["NUMEXPR_MAX_THREADS"] = str(snakemake.threads)

import lightning as L
import pst
import torch
from lightning.pytorch.callbacks import EarlyStopping, ModelCheckpoint
from torch.utils.data import DataLoader
from torchmetrics.functional.classification import multiclass_accuracy
from tqdm import tqdm

from CheckAMG.scripts.denovo.constants import NUM_CLASSES
from CheckAMG.scripts.denovo.data import GenomeDataModule
from CheckAMG.scripts.denovo.model import CheckAMGPST
from CheckAMG.scripts.denovo.utils import _resolve_accelerator, load_model, find_best_checkpoint

from CheckAMG.scripts.common.snakemake_logging import init_snakemake_logger, emit_step_banner
from CheckAMG.scripts.common.mem_limit import set_memory_limit

import pprint

log_file = snakemake.params.log
logger = init_snakemake_logger(snakemake, redirect_streams=False)
emit_step_banner(title=f"Step {snakemake.params.step[0]}/{snakemake.params.step[1]}: Train a CheckAMG-PST model", log_path=log_file, append_raw=False)
set_memory_limit(snakemake.resources.mem)
logger.debug(f"Memory limit set to {snakemake.resources.mem} GB")

import warnings
from tables import UnclosedFileWarning
warnings.filterwarnings("ignore", category=UnclosedFileWarning, module="tables")
warnings.filterwarnings("ignore", category=UserWarning, module="lightning")
warnings.filterwarnings("ignore", message=".*LeafSpec.*")

@torch.no_grad()
def test_loop(
    model: CheckAMGPST, test_loader: DataLoader, outdir: Path, threshold: float = 0.5
):
    model.eval()

    # used to decompose the multiclass label into a binary 2-label setup
    # col 0 = AVG, col 1 = viral
    decomp = torch.tensor(
        [
            [False, False],
            [False, True],
            [True, False],
            [True, True],
        ],
        device=model.device,
    )

    y_preds: list[torch.Tensor] = []
    y_trues: list[torch.Tensor] = []

    batch: pst.GenomeGraphBatch
    for batch in tqdm(
        test_loader,
        desc="Evaluating test set",
        unit="batch",
        file=stdout,
    ):
        batch = batch.to(model.device)  # type: ignore

        # shape: [n_proteins, 4 classes]
        y_pred = model.predict_step(batch)
        y_preds.append(y_pred)
        y_trues.append(batch.y)  # type: ignore

    y_pred = torch.cat(y_preds)
    y_true = torch.cat(y_trues)

    overall_accuracy = multiclass_accuracy(
        preds=y_pred,
        target=y_true,
        num_classes=NUM_CLASSES,
        average="macro",
    )

    y_pred_decomp = decomp[y_pred]
    y_true_decomp = decomp[y_true]

    y_pred_avg = y_pred_decomp[:, 0]
    y_true_avg = y_true_decomp[:, 0]

    y_pred_viral = y_pred_decomp[:, 1]
    y_true_viral = y_true_decomp[:, 1]

    avg_accuracy = multiclass_accuracy(
        preds=y_pred_avg,
        target=y_true_avg,
        num_classes=2,
        average="macro",
    )

    viral_accuracy = multiclass_accuracy(
        preds=y_pred_viral,
        target=y_true_viral,
        num_classes=2,
        average="macro",
    )

    with outdir.joinpath("test_predictions.tsv").open("w") as fp:
        header = [
            "protein_id",
            "protein_is_AVG_PRED",
            "protein_is_AVG_TRUE",
            "protein_is_viral_PRED",
            "protein_is_viral_TRUE",
        ]
        header_line = "\t".join(header)
        fp.write(f"{header_line}\n")

        for i, values in enumerate(
            zip(y_pred_avg, y_true_avg, y_pred_viral, y_true_viral)
        ):
            avg_pred, avg_true, viral_pred, viral_true = map(
                lambda x: x.item(), values
            )
            fp.write(
                f"{i}\t{avg_pred:.4f}\t{avg_true}\t{viral_pred:.4f}\t{viral_true}\n"
            )

    with outdir.joinpath("test_summary.tsv").open("w") as fp:
        fp.write("Label\tBalanced Accuracy\n")
        fp.write(f"protein_is_AVG\t{avg_accuracy:.4f}\n")
        fp.write(f"protein_is_viral\t{viral_accuracy:.4f}\n")
        fp.write(f"overall\t{overall_accuracy:.4f}\n")

    stdout.write(
        f"Test set AVG classification balanced accuracy: {avg_accuracy:.4f}\n"
    )
    stdout.write(
        f"Test set Viral gene classification balanced accuracy: {viral_accuracy:.4f}\n"
    )
    stdout.write(
        f"Overall Test set classification balanced accuracy: {overall_accuracy:.4f}\n"
    )


def setup_model(
        model_ckpt, lr,
        margin, max_ap_pairs,
        positive_mining_strategy,
        negative_mining_strategy,
        opposite_class_negatives,
        knn, context_size,
        dropout,
        verbose
        ) -> CheckAMGPST:
    model = CheckAMGPST.from_pretrained(
        model_ckpt,
        optimizer={"lr": lr},
        layer_dropout=dropout,
        loss={
            "margin": margin,
            "max_ap_pairs": max_ap_pairs,
            "positive_mining_strategy": positive_mining_strategy,
            "negative_mining_strategy": negative_mining_strategy,
            "opposite_class_negatives": opposite_class_negatives,
            "context_size": context_size,
        },
        knn=knn,
        verbose=verbose,
    )

    return model


def setup_datamodule(train_file, test_file, class_weighting) -> GenomeDataModule:
    datamodule = GenomeDataModule(
        train_file=train_file,
        test_file=test_file,
        class_weighting=class_weighting,  # type: ignore
    )

    datamodule.setup()

    return datamodule


def train(
        model: CheckAMGPST,
        datamodule: GenomeDataModule,
        batch_size: int,
        accelerator: str,
        max_epochs: int,
        logdir: str | Path,
        devices: int = 1,
        ) -> L.Trainer:
    logger.debug(f"Specified batch size: {batch_size}")
    train_loader = datamodule.train_dataloader(batch_size=batch_size)
    val_loader = datamodule.val_dataloader(batch_size=batch_size)
    logger.debug(f"Training dataloader has {len(train_loader)} batches.")
    logger.debug(f"Validation dataloader has {len(val_loader)} batches.")

    logger.debug(f"Positive mining: {model.criterion[0].positive_mining_strat}")  # type: ignore
    logger.debug(f"Negative mining: {model.criterion[0].negative_mining_strat}")  # type: ignore
    logger.debug(f"Opposite-class negatives: {model.criterion[0].opposite_class_negatives}")  # type: ignore

    monitor_value = "val_loss"
    callbacks = [
        ModelCheckpoint(
            filename="{epoch}_{train_loss:.3f}",
            save_last=True,
            save_top_k=1,
            every_n_epochs=1,
            monitor=monitor_value,
            mode="min",
        ),
        EarlyStopping(
            monitor=monitor_value,
            patience=5,
            verbose=True,
            mode="min",
            check_finite=True,
            strict=True,
            stopping_threshold=0.1,
            min_delta=0.01,
        ),
        # BatchMonitor(logdir),
    ]

    _detected_accelerator = _resolve_accelerator(accelerator)
    if _detected_accelerator == "cpu":
        torch.set_num_threads(devices)
        devices = 1

    trainer = L.Trainer(
        max_epochs=max_epochs,
        precision="bf16-mixed",
        default_root_dir=logdir,
        callbacks=callbacks,
        accelerator=accelerator,  # pytorch lightning does something more sophisticated for auto
        devices=devices,
    )
    trainer.fit(model, train_dataloaders=train_loader, val_dataloaders=val_loader)

    for stage, skip_summary in model._skipped_batches.items():
        if skip_summary:
            logger.info(
                f"Number of batches skipped during {stage=} at each epoch: {dict(skip_summary)}"
            )

    return trainer


def load_best_model(trainer: L.Trainer, accelerator: str) -> CheckAMGPST:
    if trainer.log_dir is None:
        raise ValueError("Trainer has no checkpointing directory.")

    ckpt = find_best_checkpoint(trainer.log_dir)
    return load_model(ckpt, device=accelerator)


def main():
    # Inputs
    train_file = snakemake.params.train_file
    test_file = snakemake.params.test_file
    model_ckpt = snakemake.params.model_ckpt

    # Outputs
    logdir = snakemake.params.logdir

    # Parameters
    batch_size = snakemake.params.batch_size
    lr = snakemake.params.learning_rate
    max_epochs = snakemake.params.max_epochs
    dropout = snakemake.params.dropout
    margin = snakemake.params.margin
    knn = snakemake.params.knn
    max_ap_pairs = snakemake.params.max_ap_pairs
    positive_mining_strategy = snakemake.params.positive_mining_strategy
    negative_mining_strategy = snakemake.params.negative_mining_strategy
    opposite_class_negatives = snakemake.params.opposite_class_negatives
    class_weighting = snakemake.params.class_weighting
    context_size = snakemake.params.context_size

    # Resources
    accelerator = snakemake.params.accelerator
    devices = snakemake.params.devices
    mem = snakemake.resources.mem

    # Other
    verbose = snakemake.params.debug
    seed = snakemake.params.seed

    L.seed_everything(seed)

    model = setup_model(
        model_ckpt, lr,
        margin, max_ap_pairs,
        positive_mining_strategy,
        negative_mining_strategy,
        opposite_class_negatives,
        knn, context_size,
        dropout,
        verbose
        )

    logger.info(f"Model initialized with {sum(p.numel() for p in model.parameters())} parameters.")
    logger.info(f"Model config:\n{pprint.pformat(model.config.to_dict())}")
    datamodule = setup_datamodule(train_file, test_file, class_weighting)

    trainer = train(
        model, datamodule,
        batch_size, accelerator, max_epochs,
        logdir, devices
        )

    _trainer_log_dir = trainer.log_dir

    datamodule.close("train")

    if datamodule.test_dataset is not None:
        # reload model from the best checkpoint for test set evaluation
        model = load_best_model(trainer, accelerator)

        logger.warning(
            "Doing an approximate label propagation evaluation for the input test dataset. A more holistic approach is to the use the inference workflow for model evaluation instead."
        )
        outdir = Path(_trainer_log_dir or logdir)

        if not outdir.exists():
            outdir.mkdir(parents=True)

        test_loader = datamodule.test_dataloader(batch_size=batch_size)

        test_loop(
            model=model,
            test_loader=test_loader,
            outdir=outdir,
            threshold=0.5,
        )

        datamodule.close("test")


if __name__ == "__main__":
    main()
