# CheckAMG/scripts/common/set_up_snakemake.py

import errno
import os
import platform
import select
import subprocess
import sys
from importlib.resources import files as resource_files
import logging
from CheckAMG.scripts.common.checkAMG_ASCII import ASCII
from CheckAMG.scripts.common.snakemake_logging import emit_step_banner
from CheckAMG.scripts.common.runner_logging import append_raw_to_log


def _parent_stdout_is_tty() -> bool:
    # Check the *real* underlying stdout, not the TeeStream wrapper.
    try:
        return bool(sys.__stdout__.isatty())
    except Exception:
        return False


def _stream_via_pty(snakemake_command):
    """Run snakemake under a pseudo-TTY so tqdm/rich in child scripts
    detect a terminal and emit \\r-based in-place updates. The TeeStream
    on sys.stdout passes \\r through to the screen and rewrites it as \\n
    in the log file.
    """
    import pty

    master_fd, slave_fd = pty.openpty()
    try:
        process = subprocess.Popen(
            snakemake_command,
            stdout=slave_fd,
            stderr=slave_fd,
            close_fds=True,
        )
    except Exception:
        os.close(master_fd)
        os.close(slave_fd)
        raise
    # Close the slave end in the parent; only the child writes to it.
    os.close(slave_fd)

    try:
        while True:
            ready, _, _ = select.select([master_fd], [], [], 0.2)
            if ready:
                try:
                    chunk = os.read(master_fd, 4096)
                except OSError as e:
                    # PTY EIO is the expected signal that the child closed
                    # its end of the tty (process exited).
                    if e.errno in (errno.EIO,):
                        break
                    raise
                if not chunk:
                    break
                sys.stdout.write(chunk.decode("utf-8", errors="replace"))
                sys.stdout.flush()
            elif process.poll() is not None:
                # Process exited; drain anything still buffered.
                try:
                    while True:
                        chunk = os.read(master_fd, 4096)
                        if not chunk:
                            break
                        sys.stdout.write(chunk.decode("utf-8", errors="replace"))
                        sys.stdout.flush()
                except OSError:
                    pass
                break
    finally:
        try:
            os.close(master_fd)
        except OSError:
            pass

    process.wait()
    return process.returncode


def _stream_via_pipe(snakemake_command):
    """Fallback for non-TTY parents: line-buffered pipe (the legacy path).
    Each tqdm update arrives as its own line, which is fine for files/CI.
    """
    process = subprocess.Popen(
        snakemake_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    if process.stdout:
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
    process.wait()
    return process.returncode

def get_paths(module_dir_name):
    """Access the packaged files and scripts directories"""
    try:
        files_dir = str(resource_files("CheckAMG").joinpath("files"))
    except ModuleNotFoundError as e:
        raise RuntimeError("Package data not found. Is 'CheckAMG/files' included in your package?") from e
    try:
        steps_dir = str(resource_files("CheckAMG").joinpath("scripts", module_dir_name))
    except ModuleNotFoundError as e:
        raise RuntimeError(f"Package data not found. Is 'CheckAMG/scripts/{module_dir_name}' included in your package?") from e
    try:
        workflow_path = str(resource_files("CheckAMG").joinpath("scripts", module_dir_name, "workflow.smk"))
    except ModuleNotFoundError as e:
        raise RuntimeError(f"Snakemake workflow not found. Is 'CheckAMG/scripts/{module_dir_name}/workflow.smk' included in your package?") from e
    try:
        common_dir = str(resource_files("CheckAMG").joinpath("scripts", "common"))
    except ModuleNotFoundError as e:
        raise RuntimeError("Package data not found. Is 'CheckAMG/scripts/common' included in your package?") from e

    return files_dir, steps_dir, common_dir, workflow_path
        

def create_output_dir(output_dir):
    """Create the output directory if it doesn't already exist."""
    os.makedirs(output_dir, exist_ok=True)


def run_snakemake(module_dir_name, module_arg_name, config_path, scripts_dir, args, checkamg_version, logger, log_file_path, command, post_command_warning=None):
    """Run the Snakemake pipeline using the generated config file."""
    append_raw_to_log(log_file_path, f"{ASCII}\n")

    logger.info(f"CheckAMG version {checkamg_version}")
    logger.info(f"Starting CheckAMG {module_arg_name}...")

    current_os = platform.system()
    if current_os == "Darwin":
        logger.warning(
            f"The detected OS is {current_os}, which means no hard memory limit can be set. "
            "This should be fine, but there may be problems/crashes if you are working with very large inputs that exceed your available memory."
        )
    elif current_os == "Windows":
        logger.error(f"The detected OS is {current_os}, which is not supported. Exiting...")
        raise OSError("Windows is not supported for CheckAMG.")

    # Resolve the packaged workflow path for this module (NOT __main__.py directory)
    files_dir, steps_dir, common_dir, workflow_path = get_paths(module_dir_name)

    logger.info("Command executed: %s", command)
    if post_command_warning:
        logger.warning(post_command_warning)
    try:
        input_type = args.input_type
    except AttributeError:
        try:
            if args.query_proteins and os.path.isfile(args.query_proteins):
                input_type = "faa"
            elif args.query_bin_proteins and os.path.isfile(args.query_bin_proteins):
                input_type = "faa"
            elif args.query_contigs and os.path.isfile(args.query_contigs):
                input_type = "fna"
            elif args.query_bins and os.path.isfile(args.query_bins):
                input_type = "fna"
            elif args.query_file_esm2 and os.path.isfile(args.query_file_esm2):
                input_type = "esm2"
            elif args.query_file and os.path.isfile(args.query_file):
                input_type = "pst_embed"
            elif (args.annotate_dir and os.path.isdir(args.annotate_dir)) & (args.denovo_dir and os.path.isdir(args.denovo_dir)):
                input_type = "aggregate"
        except AttributeError:
            input_type = "unknown"
    logger.debug("The input type is %s", input_type)
    logger.debug("Using workflow: %s", workflow_path)

    devices = args.devices if hasattr(args, "devices") else args.threads
    try:
        snakemake_command = [
            "snakemake",
            "--snakefile", workflow_path,
            "--nolock",
            "--configfile", config_path,
            "--directory", args.output,
            "--cores", str(devices),
            "--rerun-triggers", "input",
            "--keep-incomplete",
            "--ignore-incomplete",
        ]

        if hasattr(args, "verbose") and args.verbose:
            snakemake_command += ["--verbose", "all"]
        elif hasattr(args, "debug") and args.debug:
            snakemake_command += ["--verbose", "all"]
        else:
            snakemake_command += ["--quiet", "all"]

        # Route snakemake's output through a pseudo-TTY when checkamg is
        # being run interactively, so child scripts using tqdm/rich emit
        # \r-based in-place progress updates that the TeeStream forwards
        # to the screen verbatim and rewrites as \n in the log file. When
        # the parent isn't on a TTY (CI, output redirected to a file),
        # fall back to a line-buffered pipe.
        if _parent_stdout_is_tty():
            returncode = _stream_via_pty(snakemake_command)
        else:
            returncode = _stream_via_pipe(snakemake_command)

        if returncode != 0:
            raise subprocess.CalledProcessError(returncode, snakemake_command)

        emit_step_banner(title=f"The CheckAMG {module_arg_name} pipeline is complete.", log_path=log_file_path, append_raw=False)

    except subprocess.CalledProcessError:
        logger.error(f"CheckAMG {module_arg_name} ended prematurely with an error!")
        raise