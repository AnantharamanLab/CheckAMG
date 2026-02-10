from pathlib import Path
from sys import stdout

import lightning as L
import pst
import torch
from lightning.pytorch.callbacks import EarlyStopping, ModelCheckpoint
from torch.utils.data import DataLoader
from torchmetrics.functional.classification import multiclass_accuracy
from tqdm import tqdm

from CheckAMG.scripts.pst_checkamg.cli import Args
from CheckAMG.scripts.pst_checkamg.constants import NUM_CLASSES
from CheckAMG.scripts.pst_checkamg.data import GenomeDataModule
from CheckAMG.scripts.pst_checkamg.embed import embed
from CheckAMG.scripts.pst_checkamg.model import CheckAMGPST
from CheckAMG.scripts.pst_checkamg.utils import _resolve_accelerator, load_model


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


def setup_model(args: Args) -> CheckAMGPST:
    model = CheckAMGPST.from_pretrained(
        args.model_ckpt,
        optimizer={"lr": args.lr},
        loss={
            "margin": args.margin,
            "max_ap_pairs": args.max_ap_pairs,
            "positive_mining_strategy": args.positive_mining_strategy,
            "negative_mining_strategy": args.negative_mining_strategy,
        },
        knn=args.knn,
        n_attn_layers=args.n_attn_layers,
        n_attn_heads=args.n_attn_heads,
        verbose=args.verbose,
    )

    return model


def setup_datamodule(args: Args) -> GenomeDataModule:
    datamodule = GenomeDataModule(
        train_file=args.train_file,
        test_file=args.test_file,
        class_weighting=args.class_weighting,  # type: ignore
    )

    datamodule.setup()

    return datamodule


def train(model: CheckAMGPST, datamodule: GenomeDataModule, args: Args) -> L.Trainer:
    train_loader = datamodule.train_dataloader(batch_size=args.batch_size)
    val_loader = datamodule.val_dataloader(batch_size=args.batch_size)

    print(f"Positive mining: {model.criterion[0].positive_mining_strat}")  # type: ignore
    print(f"Negative mining: {model.criterion[0].negative_mining_strat}")  # type: ignore

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
        # BatchMonitor(args.logdir),
    ]

    _detected_accelerator = _resolve_accelerator(args.accelerator)
    if _detected_accelerator == "cpu":
        torch.set_num_threads(args.devices)
        args.devices = 1

    trainer = L.Trainer(
        max_epochs=args.max_epochs,
        precision="bf16-mixed",
        default_root_dir=args.logdir,
        callbacks=callbacks,
        accelerator=args.accelerator,  # pytorch lightning does something more sophisticated for auto
        devices=args.devices,
    )
    trainer.fit(model, train_dataloaders=train_loader, val_dataloaders=val_loader)

    for stage, skip_summary in model._skipped_batches.items():
        if skip_summary:
            print(
                f"Number of batches skipped during {stage=} at each epoch: {dict(skip_summary)}"
            )

    return trainer


def load_best_model(trainer: L.Trainer, accelerator: str) -> CheckAMGPST:
    if trainer.log_dir is not None:
        ckptdir = Path(trainer.log_dir)
        ckpt = next(ckptdir.rglob("checkpoints/epoch*.ckpt"))

        # reload model from best checkpoint
        model = load_model(ckpt, device=accelerator)

        return model

    raise ValueError("Trainer has no checkpointing directory.")


def main():
    args = Args.parse_from_cli()
    L.seed_everything(args.seed)

    model = setup_model(args)
    datamodule = setup_datamodule(args)

    trainer = train(model, datamodule, args)

    # reload model
    if datamodule.test_dataset is not None or args.train_embed_output is not None:
        model = load_best_model(trainer, args.accelerator)

    if args.train_embed_output is not None:
        print(f"Embedding training proteins -> {args.train_embed_output}")
        embed(
            model=model,
            dataset=datamodule.train_dataset,  # type: ignore
            output=args.train_embed_output,
            batch_size=args.batch_size,
            verbose=args.verbose,
        )

    datamodule.close("train")

    if datamodule.test_dataset is not None:
        print(
            "Doing an approximate label propagation evaluation for the input test dataset. A more holistic approach is to the use the inference workflow."
        )
        if trainer.logger is None:
            outdir = Path(args.logdir)
        else:
            outdir = Path(trainer.logger.log_dir or args.logdir)

        if not outdir.exists():
            outdir.mkdir(parents=True)

        test_loader = datamodule.test_dataloader(batch_size=args.batch_size)

        test_loop(
            model=model,
            test_loader=test_loader,
            outdir=outdir,
            threshold=0.5,
        )

        datamodule.close("test")


if __name__ == "__main__":
    main()
