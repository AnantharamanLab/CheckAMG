# CheckAMG/scripts/common/runner_logging.py

from __future__ import annotations

import logging
import re
import sys
import traceback
from pathlib import Path
from typing import IO


# Matches ANSI escape sequences emitted by Rich / PyTorch Lightning / tqdm:
#   ESC [ ... final-byte          (CSI sequences, e.g. \x1B[2K, \x1B[?25l, \x1B[38;2;...m)
#   ESC ] ... BEL                 (OSC sequences)
#   ESC <single-char>             (two-char escapes)
_ANSI_ESCAPE_RE = re.compile(
    r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~]|\][^\x07]*\x07)"
)


# TeeStream

class TeeStream:
    """
    Wraps a stream (stdout or stderr) so every write goes to both the
    original stream and a log file.  Drop-in replacement for sys.stdout /
    sys.stderr; satisfies the io.TextIOBase interface well enough for
    logging, argparse, click, rich, and most third-party libraries.
    """

    def __init__(
        self,
        original: IO[str],
        log_path: Path,
        mode: str = "a",
        encoding: str = "utf-8",
    ) -> None:
        self._orig = original
        self._file = log_path.open(mode, encoding=encoding, buffering=1)
        self.is_teestream = True

    # core write interface

    def write(self, data: str) -> int:
        # Pass the raw payload (including any \r and ANSI escapes) straight
        # to the original stream so tqdm/rich progress bars can update in
        # place with colors on a tty.
        self._orig.write(data)
        self._orig.flush()
        # For the log file, sanitize the same payload:
        #   * collapse \r\n -> \n first, so PTY's ONLCR newline translation
        #     doesn't produce a blank line between every record.
        #   * then turn any remaining standalone \r (tqdm/Rich in-place
        #     update markers) into \n so each progress update becomes one
        #     readable line in the log instead of overwriting the previous.
        #   * strip ANSI escape sequences (cursor hide/show, line erase,
        #     SGR colors) so Rich-style progress bars from PyTorch Lightning
        #     leave plain text in the log instead of \x1B[2K\x1B[38;2;...m
        #     garbage.
        if "\r" in data:
            log_data = data.replace("\r\n", "\n").replace("\r", "\n")
        else:
            log_data = data
        if "\x1B" in log_data:
            log_data = _ANSI_ESCAPE_RE.sub("", log_data)
        self._file.write(log_data)
        self._file.flush()
        return len(data)

    def writelines(self, lines) -> None:
        for line in lines:
            self.write(line)

    def flush(self) -> None:
        self._orig.flush()
        self._file.flush()

    def close(self) -> None:
        try:
            self._file.flush()
            self._file.close()
        except Exception:
            pass

    # compatibility shims

    @property
    def encoding(self) -> str:
        return getattr(self._orig, "encoding", None) or "utf-8"

    @property
    def errors(self) -> str | None:
        return getattr(self._orig, "errors", "replace")

    def isatty(self) -> bool:
        # Proxy to the original stream so tqdm/rich/click decide their
        # output style based on the *real* terminal, not the tee wrapper.
        try:
            return bool(self._orig.isatty())
        except Exception:
            return False

    def fileno(self) -> int:
        # Some tools call fileno() to check whether they can use low-level
        # fd tricks. Fall back to the original fd so they don't crash
        return self._orig.fileno()

    @property
    def closed(self) -> bool:
        return self._file.closed

    def readable(self) -> bool:
        return False

    def writable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return False

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


# stream / excepthook installation

def redirect_streams_to_log(
    log_file_path: str | Path,
    *,
    capture_stdout: bool = True,
    capture_stderr: bool = True,
    redirect_logging: bool = True,
) -> None:
    """
    Replace sys.stdout and/or sys.stderr with TeeStreams that mirror every
    write to *log_file_path*.  Also installs sys.excepthook so unhandled
    exceptions are written to the log before the process exits.
        setup_runner_logger(log_path, debug) # structured logger
        redirect_streams_to_log(log_path)    # capture everything else
    """
    # If sys.stdout is already a TeeStream, don't wrap it again.
    if getattr(sys.stdout, "is_teestream", False):
        return 

    p = Path(log_file_path)
    p.parent.mkdir(parents=True, exist_ok=True)

    if capture_stdout:
        sys.stdout = TeeStream(sys.__stdout__, p)

    if capture_stderr:
        sys.stderr = TeeStream(sys.__stderr__, p)

    # Hook into the root logger so that logging calls from any logger are captured by the TeeStreams
    if redirect_logging:
        root_logger = logging.getLogger()
        # Prevent duplicate handlers if called multiple times
        if not any(isinstance(h, logging.FileHandler) and str(h.baseFilename) == str(p.absolute()) for h in root_logger.handlers):
            fmt = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
            fh = logging.FileHandler(p, mode="a")
            fh.setFormatter(fmt)
            root_logger.addHandler(fh)

    # Unhandled exceptions are printed to stderr by the default excepthook,
    # which is now a TeeStream, so they land in the log automatically.
    # Override anyway to add a clear prefix and guarantee flushing.
    _original_excepthook = sys.excepthook

    def _excepthook(exc_type, exc_value, exc_tb):
        header = "\n[UNHANDLED EXCEPTION]\n"
        tb_text = "".join(traceback.format_exception(exc_type, exc_value, exc_tb))
        sys.stderr.write(header + tb_text)
        sys.stderr.flush()

    sys.excepthook = _excepthook


def restore_streams() -> None:
    """
    Restore sys.stdout / sys.stderr to the original system streams and close
    the underlying log file handles.  Call in a finally block if you need
    clean teardown (e.g. in tests).
    """
    for attr, orig in (("stdout", sys.__stdout__), ("stderr", sys.__stderr__)):
        current = getattr(sys, attr)
        if isinstance(current, TeeStream):
            current.close()
        setattr(sys, attr, orig)


def setup_runner_logger(
    log_file_path: str | Path | None,
    debug: bool,
) -> logging.Logger:
    # Immediately install redirection so early logs are captured
    if log_file_path:
        redirect_streams_to_log(log_file_path)

    logger = logging.getLogger("CheckAMG.runner")
    logger.setLevel(logging.DEBUG if debug else logging.INFO)
    logger.propagate = False

    for h in list(logger.handlers):
        logger.removeHandler(h)

    fmt = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Only use StreamHandler. TeeStream handles the file writing.
    sh = logging.StreamHandler(sys.stdout) 
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    return logger


def append_raw_to_log(log_file_path: str | Path, text: str) -> None:
    p = Path(log_file_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("a", encoding="utf-8") as f:
        f.write(text)
        if not text.endswith("\n"):
            f.write("\n")


def print_and_append(log_file_path: str | Path, text: str) -> None:
    print(text, flush=True)
    append_raw_to_log(log_file_path, text)