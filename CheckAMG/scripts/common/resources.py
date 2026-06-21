"""Resource-allocation helpers shared by CheckAMG subcommands."""

from __future__ import annotations


def _gpu_will_be_used(accelerator: str) -> bool:
    """Return True if the effective device after Lightning's auto-resolution will be GPU."""
    if accelerator == "gpu":
        return True
    if accelerator == "auto":
        try:
            import torch
            return bool(torch.cuda.is_available())
        except Exception:
            return False
    return False


def clamp_gpu_devices(accelerator: str, devices: int) -> tuple[int, str | None]:
    """If GPU will be used and ``devices > 1``, clamp to 1 and return a warning
    message that the caller can emit at the appropriate point in the log stream.

    Returns ``(effective_devices, warning_or_None)``.

    Multi-GPU is not supported in any CheckAMG denovo/train step: the ESM2
    embedding writer is not DDP-aware (the per-rank writes race on a single
    h5 file), and the PST embedding/inference paths run only on the single
    device the model was loaded onto. Multi-GPU training would also break
    bit-identical reproducibility (bf16 all-reduce non-associativity, and
    triplet-loss mining quality depends on full-batch composition), so it is
    unsupported here.
    """
    if devices > 1 and _gpu_will_be_used(accelerator):
        return 1, (
            f"Multi-GPU is not supported (--accelerator {accelerator} --devices {devices}); "
            f"clamping to --devices 1."
        )
    return devices, None
