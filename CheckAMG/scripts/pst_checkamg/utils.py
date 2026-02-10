from pathlib import Path

import torch

from CheckAMG.scripts.pst_checkamg.model import CheckAMGPST

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
