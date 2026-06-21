# CheckAMG/scripts/common/mem_limit.py

from __future__ import annotations
import resource
from CheckAMG.scripts.common.snakemake_logging import init_snakemake_logger

def set_memory_limit(limit_gb: float | int | None) -> None:
    """
    Set a per-process address-space limit (RLIMIT_AS) in GB.
    No-op if limit_gb is None, 0, or negative.

    Note: RLIMIT_AS is best-effort and may behave differently across systems.
    """
    if limit_gb is None:
        return
    try:
        limit = float(limit_gb)
    except Exception:
        return
    if limit <= 0:
        return

    limit_bytes = int(limit * (1024 ** 3))
    try:
        soft, hard = resource.getrlimit(resource.RLIMIT_AS)
        new_soft = min(limit_bytes, hard) if hard not in (-1, resource.RLIM_INFINITY) else limit_bytes
        new_hard = hard
        resource.setrlimit(resource.RLIMIT_AS, (new_soft, new_hard))
    except Exception:
        raise RuntimeError(f"Failed to set memory limit to {limit_gb} GB (RLIMIT_AS={limit_bytes} bytes).")