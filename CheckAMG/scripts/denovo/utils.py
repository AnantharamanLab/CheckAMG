from pathlib import Path

import torch

from CheckAMG.scripts.denovo.model import CheckAMGPST
from CheckAMG.scripts.denovo.data import GenomeDataModule

CUDA_AVAILABLE = torch.cuda.is_available()


# quick check -> lightning does something more sophisticated
def _resolve_accelerator(accelerator: str) -> str:
    if accelerator == "auto":
        if CUDA_AVAILABLE:
            return "gpu"
        return "cpu"
    return accelerator


def load_model(
    ckpt: str | Path, device: torch.device | None | str = None
) -> CheckAMGPST:
    model = CheckAMGPST.from_pretrained(ckpt)

    device_is_str = isinstance(device, str)
    if device is None or (device_is_str and device == "auto"):
        device = torch.device("cuda" if CUDA_AVAILABLE else "cpu")
    elif device_is_str:
        if device == "gpu" and CUDA_AVAILABLE:
            device = torch.device("cuda")
        else:
            device = torch.device("cpu")

    model = model.to(device)

    return model.eval()


def find_best_checkpoint(logdir: str | Path) -> Path:
    logdir = Path(logdir)
    checkpoints = list(logdir.rglob("checkpoints/epoch*.ckpt"))
    if not checkpoints:
        raise FileNotFoundError(
            f"No trained model checkpoint (checkpoints/epoch*.ckpt) found under {logdir}."
        )
    # The most recently written checkpoint corresponds to the latest training run.
    return max(checkpoints, key=lambda p: p.stat().st_mtime)


def _extract_label_from_datamodule(datamodule: GenomeDataModule) -> torch.Tensor:
    ptr: torch.Tensor = datamodule.train_dataset.dataset.scaffold_ptr  # type: ignore
    full_y: torch.Tensor = datamodule.train_dataset.dataset.y  # type: ignore
    scaffold_indices = datamodule.train_dataset.indices  # type: ignore

    train_y: list[torch.Tensor] = []
    for idx in scaffold_indices:
        train_y.append(full_y[ptr[idx] : ptr[idx + 1]])

    return torch.cat(train_y)