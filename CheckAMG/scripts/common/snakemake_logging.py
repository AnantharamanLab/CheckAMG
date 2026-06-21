# CheckAMG/scripts/common/snakemake_logging.py

from __future__ import annotations

import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

@dataclass(frozen=True)
class LoggerConfig:
    name: str = "CheckAMG"
    level: int = logging.INFO
    log_path: Optional[Path] = None
    also_stdout: bool = True
    mode: str = "a"
    fmt: str = "%(asctime)s | %(levelname)s | %(message)s"
    datefmt: str = "%Y-%m-%d %H:%M:%S"

def _get_bool(obj: Any, key: str, default: bool = False) -> bool:
    try:
        v = getattr(obj, key)
        if isinstance(v, bool):
            return v
        if v is None:
            return default
        if isinstance(v, (int, float)):
            return bool(v)
        if isinstance(v, str):
            return v.strip().lower() in {"1", "true", "t", "yes", "y", "on"}
        return default
    except Exception:
        return default

def _as_path(p: Any) -> Optional[Path]:
    if p is None:
        return None
    try:
        s = str(p).strip()
        if not s:
            return None
        return Path(s)
    except Exception:
        return None

def get_snakemake_log_path(
    snakemake: Any,
    *,
    log_key: str = "log",
) -> Optional[Path]:

    try:
        params = getattr(snakemake, "params", None)
        if params is not None and hasattr(params, log_key):
            p = _as_path(getattr(params, log_key))
            if p is not None:
                return p
    except Exception:
        pass

    try:
        snk_log = getattr(snakemake, "log", None)
        if isinstance(snk_log, (list, tuple)) and snk_log:
            return _as_path(snk_log[0])
        return _as_path(snk_log)
    except Exception:
        return None

def init_snakemake_logger(
    snakemake: Any,
    name: str = "CheckAMG",
    *,
    log_key: str = "log",
    debug_param: str = "debug",
    level: Optional[int] = None,
    also_stdout: bool = True,
    mode: str = "a",
    fmt: str = "%(asctime)s | %(levelname)s | %(message)s",
    datefmt: str = "%Y-%m-%d %H:%M:%S",
    redirect_streams: bool = True,
) -> logging.Logger:
    debug = False
    try:
        params = getattr(snakemake, "params", None)
        if params is not None:
            debug = _get_bool(params, debug_param, False)
    except Exception:
        debug = False

    eff_level = level if level is not None else (logging.DEBUG if debug else logging.INFO)
    log_path = get_snakemake_log_path(snakemake, log_key=log_key)

    cfg = LoggerConfig(
        name=name,
        level=eff_level,
        log_path=log_path,
        also_stdout=also_stdout,
        mode=mode,
        fmt=fmt,
        datefmt=datefmt,
    )
    logger = _build_logger(cfg)

    # Mirror all raw stdout/stderr (tqdm, rich, third-party libs) to the log file
    # Must come AFTER _build_logger so the FileHandler is attached first
    if redirect_streams and log_path is not None:
        from CheckAMG.scripts.common.runner_logging import redirect_streams_to_log
        redirect_streams_to_log(log_path)

    return logger

def _build_logger(cfg: LoggerConfig) -> logging.Logger:
    logger = logging.getLogger(cfg.name)
    logger.setLevel(cfg.level)
    logger.propagate = False

    for h in list(logger.handlers):
        logger.removeHandler(h)

    formatter = logging.Formatter(cfg.fmt, datefmt=cfg.datefmt)

    if cfg.also_stdout:
        sh = logging.StreamHandler(sys.stdout)
        sh.setLevel(cfg.level)
        sh.setFormatter(formatter)
        logger.addHandler(sh)

    return logger

def format_step_banner(
    title: str,
    *,
    width: int = 80,
    border: str = "=",
    pad: str = " ",
) -> str:
    if width < 20:
        width = 20
    t = str(title).strip()
    top = border * width
    mid = t.center(width, pad)
    bot = border * width
    return "\n".join([top, mid, bot])

def append_banner_to_file(log_path: Path, banner: str) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as f:
        f.write(banner + "\n")

def emit_step_banner(
    title: str,
    *,
    snakemake: Any | None = None,
    log_path: str | Path | None = None,
    width: int = 80,
    border: str = "=",
    stdout: bool = True,
    append_raw: bool = True,
    log_key: str = "log",
) -> None:
    """
    Emit a fixed-width banner without logger prefixes.

    - If stdout=True, prints banner to stdout (no timestamps/levels).
    - If append_raw=True, appends banner to the log file (no timestamps/levels).

    Pass either snakemake=snakemake (to resolve log path) or log_path directly.
    """
    banner = format_step_banner(title, width=width, border=border)

    if stdout:
        print(banner, file=sys.stdout, flush=True)

    # Use the attribute check instead of isinstance
    if getattr(sys.stdout, "is_teestream", False):
        # The TeeStream has already written the banner to the file via 'print'
        return 

    if not append_raw:
        return

    lp: Optional[Path] = None
    if log_path is not None:
        lp = Path(log_path)
    elif snakemake is not None:
        lp = get_snakemake_log_path(snakemake, log_key=log_key)

    if lp is None:
        return

    try:
        append_banner_to_file(lp, banner)
    except Exception:
        return