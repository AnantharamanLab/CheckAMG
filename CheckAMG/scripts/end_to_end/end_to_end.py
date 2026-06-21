#!/usr/bin/env python3

import logging
import os
import re
import sys
from pathlib import Path

from CheckAMG.scripts.common.checkAMG_ASCII import ASCII
from CheckAMG.scripts.common.runner_logging import (
    _ANSI_ESCAPE_RE,
    append_raw_to_log,
    restore_streams,
)
from CheckAMG.scripts.common.snakemake_logging import format_step_banner
from CheckAMG.scripts.common.set_up_snakemake import create_output_dir, run_snakemake
from CheckAMG.scripts.annotate import annotate
from CheckAMG.scripts.denovo import denovo
from CheckAMG.scripts.aggregate import aggregate

_STEP_RE = re.compile(r"Step\s+(\d+)\s*/\s*(\d+)\s*:\s*(.+?)\s*$")


class _MonitorController:
    """Holds the currently active module log file and the step callback.

    The streams below mirror everything that the running module would normally
    print to its per-module log file (preserving the standalone module logs),
    while the orchestrator drives a single high-level end-to-end log. Detailed
    module output is intentionally not echoed to the terminal.
    """

    def __init__(self):
        self._fh = None
        self._mirror_fh = None
        self.on_step = None

    def set_mirror(self, log_path):
        """Mirror all per-module output into the end-to-end log (debug mode)."""
        p = Path(log_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        self._mirror_fh = p.open("a", encoding="utf-8", buffering=1)

    def set_module_log(self, log_path):
        if self._fh is not None:
            try:
                self._fh.flush()
                self._fh.close()
            except Exception:
                pass
        p = Path(log_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        self._fh = p.open("a", encoding="utf-8", buffering=1)

    def write_log(self, data):
        if "\r" in data:
            data = data.replace("\r\n", "\n").replace("\r", "\n")
        if "\x1B" in data:
            data = _ANSI_ESCAPE_RE.sub("", data)
        for fh in (self._fh, self._mirror_fh):
            if fh is not None:
                fh.write(data)
                fh.flush()

    def scan_line(self, line):
        if self.on_step is None:
            return
        clean = _ANSI_ESCAPE_RE.sub("", line)
        m = _STEP_RE.search(clean)
        if m:
            self.on_step(int(m.group(1)), int(m.group(2)), m.group(3).strip())

    def close(self):
        for attr in ("_fh", "_mirror_fh"):
            fh = getattr(self, attr)
            if fh is not None:
                try:
                    fh.flush()
                    fh.close()
                except Exception:
                    pass
                setattr(self, attr, None)


class _MonitorStream:
    """Stdout/stderr replacement used during an end-to-end run.

    Marks itself as a TeeStream so the per-module redirect_streams_to_log calls
    short-circuit and do not re-wrap the streams or attach their own file
    handlers. Writes are routed to the active module's log file (and scanned for
    step banners on stdout), but never to the real terminal.
    """

    def __init__(self, original, controller, scan):
        self._orig = original
        self._ctrl = controller
        self._scan = scan
        self._linebuf = ""
        self.is_teestream = True

    def write(self, data):
        self._ctrl.write_log(data)
        if self._scan:
            self._linebuf += data
            parts = re.split(r"[\r\n]", self._linebuf)
            self._linebuf = parts.pop()
            for line in parts:
                self._ctrl.scan_line(line)
        return len(data)

    def writelines(self, lines):
        for line in lines:
            self.write(line)

    def flush(self):
        pass

    def close(self):
        pass

    @property
    def encoding(self):
        return getattr(self._orig, "encoding", None) or "utf-8"

    @property
    def errors(self):
        return getattr(self._orig, "errors", "replace")

    def isatty(self):
        try:
            return bool(self._orig.isatty())
        except Exception:
            return False

    def fileno(self):
        return self._orig.fileno()

    @property
    def closed(self):
        return False

    def readable(self):
        return False

    def writable(self):
        return True

    def seekable(self):
        return False


def _setup_e2e_logger(log_file_path, debug, terminal):
    logger = logging.getLogger("CheckAMG.end_to_end")
    logger.setLevel(logging.DEBUG if debug else logging.INFO)
    logger.propagate = False
    for h in list(logger.handlers):
        logger.removeHandler(h)

    fmt = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S")

    fh = logging.FileHandler(log_file_path, mode="a")
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    sh = logging.StreamHandler(terminal)
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    return logger


def _annotate_protein_inputs(annotate_output):
    """Return (single_faa, bin_faa_dir) for the proteins annotate produced and
    filtered, each only if present and non-empty.

    These are the CDS-filtered translations under filtered_faa_by_cds (the raw
    pyrodigal-gv directory is cleaned up after filtering). The same location is
    used for both nucleotide and protein inputs.
    """
    fdir = os.path.join(annotate_output, "wdir", "filtered_input", "filtered_faa_by_cds")
    single = os.path.join(fdir, "single_contig_proteins.faa")
    bins = os.path.join(fdir, "bin_proteins")

    single_ok = os.path.isfile(single) and os.path.getsize(single) > 0

    bins_ok = False
    if os.path.isdir(bins):
        for f in os.listdir(bins):
            fp = os.path.join(bins, f)
            if os.path.isfile(fp) and os.path.getsize(fp) > 0:
                bins_ok = True
                break

    return (single if single_ok else None, bins if bins_ok else None)


def _count_done(module_output, cap):
    """Number of completed snakemake step markers for a module (capped)."""
    d = os.path.join(module_output, "snakemake")
    if not os.path.isdir(d):
        return 0
    n = sum(1 for f in os.listdir(d) if f.endswith(".done"))
    return min(n, cap)


def run_end_to_end(common, annotate_ns, denovo_ns, aggregate_ns, scripts_dir, checkamg_version, command):
    output = os.path.abspath(common.output)
    create_output_dir(output)

    annotate_ns.output = os.path.join(output, "01_annotate")
    denovo_ns.output = os.path.join(output, "02_denovo")
    aggregate_ns.output = os.path.join(output, "03_aggregate")
    create_output_dir(annotate_ns.output)
    create_output_dir(denovo_ns.output)
    create_output_dir(aggregate_ns.output)

    annotate_log = os.path.join(annotate_ns.output, "CheckAMG_annotate.log")
    denovo_log = os.path.join(denovo_ns.output, "CheckAMG_denovo.log")
    aggregate_log = os.path.join(aggregate_ns.output, "CheckAMG_aggregate.log")
    e2e_log = os.path.join(output, "CheckAMG_end_to_end.log")

    annotate_steps = 11 if annotate_ns.input_type == "nucl" else 8
    denovo_steps = 4  # always fed protein FASTAs from annotate (faa query mode)
    aggregate_steps = 1
    total_steps = annotate_steps + denovo_steps + aggregate_steps

    terminal = sys.__stdout__
    controller = _MonitorController()
    monitor_out = _MonitorStream(sys.__stdout__, controller, scan=True)
    monitor_err = _MonitorStream(sys.__stderr__, controller, scan=False)
    sys.stdout = monitor_out
    sys.stderr = monitor_err

    append_raw_to_log(e2e_log, f"{ASCII}\n")
    logger = _setup_e2e_logger(e2e_log, common.debug, terminal)

    # In debug mode, mirror every module's full output (including its debug
    # messages) into the end-to-end log, not just the per-step summaries.
    if common.debug:
        controller.set_mirror(e2e_log)

    def make_on_step(module_label, module_steps):
        seen = set()

        def on_step(n, total, title):
            if n in seen:
                return
            seen.add(n)
            logger.info(f"[{module_label}] Step {n}/{module_steps}: {title}")

        return on_step

    def emit_banner(title):
        banner = format_step_banner(title)
        append_raw_to_log(e2e_log, banner)
        print(banner, file=terminal, flush=True)

    def start_module(idx, label, module_steps, module_output, module_log, extra_lines=()):
        controller.set_module_log(module_log)
        controller.on_step = make_on_step(label, module_steps)
        unit = "step" if module_steps == 1 else "steps"
        emit_banner(f"Module {idx}/3: {label} ({module_steps} {unit})")
        logger.info(f"  outputs -> {module_output}")
        logger.info(f"  module log -> {module_log}")
        for line in extra_lines:
            logger.info(line)
        done_before = _count_done(module_output, module_steps)
        if done_before >= module_steps:
            logger.info(f"  all {module_steps} {unit} already complete; skipping (resumed run).")
        elif done_before > 0:
            logger.info(
                f"  {done_before}/{module_steps} steps already complete; resuming at step {done_before + 1}."
            )

    def finish_module(label):
        logger.info(f"CheckAMG {label} is complete.")

    try:
        emit_banner("Starting CheckAMG end-to-end")
        logger.info(f"CheckAMG version {checkamg_version}")
        logger.info("Command executed: %s", command)
        logger.info(f"Output directory: {output}")
        logger.info(f"End-to-end log: {e2e_log}")
        if common.debug:
            logger.debug(f"Annotate input type: {annotate_ns.input_type}")
            logger.debug(
                f"Planned steps: annotate={annotate_steps}, de-novo={denovo_steps}, "
                f"aggregate={aggregate_steps} (total={total_steps})"
            )

        total_done = (
            _count_done(annotate_ns.output, annotate_steps)
            + _count_done(denovo_ns.output, denovo_steps)
            + _count_done(aggregate_ns.output, aggregate_steps)
        )
        if total_done > 0:
            logger.info(f"Resuming a previous run: {total_done}/{total_steps} steps already complete.")

        # 01 annotate
        start_module(1, "annotate", annotate_steps, annotate_ns.output, annotate_log)
        config_path, mod_logger, log_file_path = annotate.generate_config(annotate_ns)
        run_snakemake(
            module_dir_name="annotate",
            module_arg_name="annotate",
            config_path=config_path,
            scripts_dir=scripts_dir,
            args=annotate_ns,
            checkamg_version=checkamg_version,
            logger=mod_logger,
            log_file_path=log_file_path,
            command=command,
        )
        finish_module("annotate")

        # Hand off the proteins annotate produced and filtered to de-novo
        single_faa, bin_faa_dir = _annotate_protein_inputs(annotate_ns.output)
        denovo_ns.query_proteins = single_faa
        denovo_ns.query_bin_proteins = bin_faa_dir
        denovo_ns.query_contigs = None
        denovo_ns.query_bins = None

        if not denovo_ns.query_proteins and not denovo_ns.query_bin_proteins:
            raise FileNotFoundError(
                "No proteins available from the annotate step to run de-novo. "
                f"Check the annotate outputs and log at {annotate_log}."
            )

        # 02 de-novo
        start_module(
            2, "de-novo", denovo_steps, denovo_ns.output, denovo_log,
            extra_lines=["  using proteins produced by the annotate module as de-novo query"],
        )
        config_path, mod_logger, log_file_path, gpu_clamp_warning = denovo.generate_config(denovo_ns)
        run_snakemake(
            module_dir_name="denovo",
            module_arg_name="de-novo",
            config_path=config_path,
            scripts_dir=scripts_dir,
            args=denovo_ns,
            checkamg_version=checkamg_version,
            logger=mod_logger,
            log_file_path=log_file_path,
            command=command,
            post_command_warning=gpu_clamp_warning,
        )
        finish_module("de-novo")

        # 03 aggregate
        aggregate_ns.annotate_dir = annotate_ns.output
        aggregate_ns.denovo_dir = denovo_ns.output
        start_module(3, "aggregate", aggregate_steps, aggregate_ns.output, aggregate_log)
        config_path, mod_logger, log_file_path = aggregate.generate_config(aggregate_ns)
        run_snakemake(
            module_dir_name="aggregate",
            module_arg_name="aggregate",
            config_path=config_path,
            scripts_dir=scripts_dir,
            args=aggregate_ns,
            checkamg_version=checkamg_version,
            logger=mod_logger,
            log_file_path=log_file_path,
            command=command,
        )
        finish_module("aggregate")

        emit_banner("The CheckAMG end-to-end pipeline is complete.")
        logger.info(f"Aggregated results -> {aggregate_ns.output}")
    except Exception:
        logger.error("CheckAMG end-to-end ended prematurely with an error!")
        raise
    finally:
        controller.close()
        restore_streams()
