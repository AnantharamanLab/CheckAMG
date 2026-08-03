#!/usr/bin/env python3

import os
import shutil
import tarfile
from pathlib import Path

import requests
from packaging.version import InvalidVersion, Version
from tqdm import tqdm

from CheckAMG.scripts.common.runner_logging import setup_runner_logger
logger = setup_runner_logger(log_file_path=None, debug=False)

# DB URL selection rules:
# - Rules are evaluated top-to-bottom; first match wins.
# - Use open-ended ranges (max=None) so newer CheckAMG versions keep using the same DB unless changed.
# - "url" is the database required by 'checkamg annotate'.
# - "denovo_url" is the larger database required by 'checkamg denovo'; it is extracted into its own
#   folder under dest, alongside (not inside) the annotate DB folder. Omit when no de-novo DB exists.
# - Can optionally set tar_name / denovo_tar_name explicitly; otherwise inferred from the URL basename.
# - db_version is a stable identifier that can be an override (e.g. "checkamg_db_v1_20260128").
# - denovo_db_version is the matching stable identifier for the de-novo DB; required whenever
#   "denovo_url" is set so the de-novo DB is version-pinned alongside the annotate DB.
CHECKAMG_DB_URL_RULES = [
    {
        "min_version": "1.1",
        "max_version": None, # open-ended
        "db_version": "CheckAMG_annotate_db_v1.1_20260316",
        "url": "https://zenodo.org/records/21776005/files/CheckAMG_annotate_db_v1.1_20260316.tar.gz?download=1",
        "denovo_db_version": "CheckAMG_denovo_db_v1.1_20260714",
        "denovo_url": "https://zenodo.org/records/21776005/files/CheckAMG_denovo_db_v1.1_20260714.tar.gz?download=1",
    },
    {
        "min_version": "1.0",
        "max_version": "1.0",
        "db_version": "CheckAMG_annotate_db_v1.1_20260316",
        "url": "https://zenodo.org/records/20750360/files/CheckAMG_annotate_db_v1.1_20260316.tar.gz?download=1",
        "denovo_db_version": "checkamg_denovo_db_v1.1_20260316",
        "denovo_url": "https://zenodo.org/records/20750360/files/CheckAMG_denovo_db_v1_20260605.tar.gz?download=1",
        # "tar_name": "CheckAMG_db_v1.1_20260316.tar.gz", # optional
        # "denovo_tar_name": "CheckAMG_denovo_db_v1.1_20260316.tar.gz", # optional
    },
    {
        "min_version": "0.7.0",
        "max_version": "0.10.0", # inclusive
        "db_version": "checkamg_db_v1_20260128",
        "url": "https://zenodo.org/records/18407279/files/CheckAMG_db_v1_20260128.tar.gz?download=1",
        # "tar_name": "CheckAMG_db_v1_20260128.tar.gz", # optional
    },
]


def _parse_version(v: str) -> Version:
    try:
        return Version(v)
    except InvalidVersion as e:
        raise RuntimeError(f"Invalid version string: {v}") from e

def get_db_rule_for_db_version(db_version: str) -> dict:
    matches = [r for r in CHECKAMG_DB_URL_RULES if r.get("db_version") == db_version]
    if not matches:
        available = ", ".join(sorted({r.get("db_version", "unknown") for r in CHECKAMG_DB_URL_RULES}))
        raise RuntimeError(f"Unknown db_version '{db_version}'. Available: {available}")
    if len(matches) > 1:
        raise RuntimeError(f"Ambiguous db_version '{db_version}': multiple rules match.")
    rule = matches[0]
    if "url" not in rule:
        raise RuntimeError(f"DB rule missing 'url' for db_version '{db_version}': {rule}")
    return rule

def get_db_rule_for_checkamg(checkamg_version: str) -> dict:
    v = _parse_version(checkamg_version)

    for rule in CHECKAMG_DB_URL_RULES:
        vmin = _parse_version(rule["min_version"])
        vmax = _parse_version(rule["max_version"]) if rule.get("max_version") is not None else None

        if v < vmin:
            continue
        if vmax is not None and v > vmax:
            continue

        if "url" not in rule:
            raise RuntimeError(f"DB rule missing 'url': {rule}")
        return rule

    raise RuntimeError(f"No DB URL rule matches CheckAMG version: {checkamg_version}")

def _tar_name_from_url(url: str) -> str:
    name = url.split("?", 1)[0].rstrip("/").split("/")[-1]
    if not name or not name.endswith(".tar.gz"):
        raise RuntimeError(f"Could not infer .tar.gz filename from URL: {url}")
    return name

def _extract_dir_name_from_tar_name(tar_name: str) -> str:
    if tar_name.endswith(".tar.gz"):
        return tar_name[: -len(".tar.gz")]
    if tar_name.endswith(".tgz"):
        return tar_name[: -len(".tgz")]
    return tar_name

def _safe_extract_tar(tar: tarfile.TarFile, dest_dir: Path) -> None:
    dest_dir = dest_dir.resolve()
    for member in tar.getmembers():
        member_path = (dest_dir / member.name).resolve()
        if not str(member_path).startswith(str(dest_dir) + os.sep) and member_path != dest_dir:
            raise RuntimeError(f"Unsafe tar entry detected (path traversal): {member.name}")
    tar.extractall(path=dest_dir)

def _download_file(url: str, outpath: Path, force: bool = False, chunk_size: int = 1024 * 1024) -> None:
    outpath.parent.mkdir(parents=True, exist_ok=True)

    if outpath.exists() and not force:
        logger.info(f"Tarball already exists at {outpath}. Skipping download.")
        return

    tmp = outpath.with_suffix(outpath.suffix + ".tmp")
    if tmp.exists():
        tmp.unlink()

    logger.info(f"Downloading {url}")
    with requests.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        total = r.headers.get("Content-Length")
        total_bytes = int(total) if total is not None else None

        with open(tmp, "wb") as f, tqdm(
            total=total_bytes,
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
            desc=outpath.name,
            leave=True,
        ) as pbar:
            for chunk in r.iter_content(chunk_size=chunk_size):
                if not chunk:
                    continue
                f.write(chunk)
                pbar.update(len(chunk))

    tmp.replace(outpath)
    logger.info(f"Downloaded to {outpath}")

def _all_members_share_single_topdir(tar: tarfile.TarFile) -> str | None:
    top = None
    for m in tar.getmembers():
        if not m.name or m.name in (".", "./"):
            continue
        name = m.name.lstrip("./")
        first = name.split("/", 1)[0]
        if first in ("", ".", ".."):
            return None
        if top is None:
            top = first
        elif first != top:
            return None
    return top

def _extract_tar_gz(tar_gz_path: Path, extract_dir: Path, force: bool = False) -> None:
    if extract_dir.exists():
        if force:
            shutil.rmtree(extract_dir)
        else:
            raise RuntimeError(f"Database already exists at {extract_dir}. Use --force to overwrite.")

    extract_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Extracting {tar_gz_path}")
    with tarfile.open(tar_gz_path, mode="r:gz") as tar:
        topdir = _all_members_share_single_topdir(tar)

        if topdir is None:
            # Tar has multiple top-level entries; extract as-is into extract_dir
            _safe_extract_tar(tar, extract_dir)
            return

        # Tar has a single top-level dir; strip it so contents land directly in extract_dir
        dest_dir = extract_dir.resolve()
        for member in tar.getmembers():
            if not member.name or member.name in (".", "./"):
                continue

            name = member.name.lstrip("./")
            if name == topdir:
                continue
            if not name.startswith(topdir + "/"):
                # Shouldn't happen if topdir detection is correct; be conservative.
                raise RuntimeError(f"Unexpected tar member outside topdir '{topdir}': {member.name}")

            stripped = name[len(topdir) + 1 :]
            if not stripped:
                continue

            member.name = stripped  # tarfile uses this for extraction path

            # Repeat your traversal safety check with the stripped path
            member_path = (dest_dir / member.name).resolve()
            if not str(member_path).startswith(str(dest_dir) + os.sep) and member_path != dest_dir:
                raise RuntimeError(f"Unsafe tar entry detected (path traversal): {member.name}")

            tar.extract(member, path=dest_dir)

def download_annotate_db(rule: dict, dest_dir: Path, checkamg_version: str, force: bool = False) -> Path:
    url = rule["url"]
    chosen_db_version = rule.get("db_version", "unknown")

    tar_name = rule.get("tar_name") or _tar_name_from_url(url)
    tar_path = dest_dir / tar_name

    extract_dir_name = rule.get("extract_dir") or _extract_dir_name_from_tar_name(tar_name)
    extract_dir = dest_dir / extract_dir_name

    logger.info(f"Using database version: {chosen_db_version} for CheckAMG version {checkamg_version}")
    logger.info(f"Annotate database source URL: {url}")

    if extract_dir.exists() and not force:
        logger.error(f"Database {chosen_db_version} already exists at {extract_dir}. Use --force to overwrite.")
        raise RuntimeError(f"Database {chosen_db_version} already exists at {extract_dir}. Use --force to overwrite.")

    _download_file(url, tar_path, force=force)
    _extract_tar_gz(tar_path, extract_dir, force=force)

    if tar_path.exists():
        tar_path.unlink()

    logger.info(f"Annotate database download and extraction complete. Location: {extract_dir}")
    return extract_dir

def download_denovo_db(rule: dict, dest_dir: Path, force: bool = False) -> Path | None:
    denovo_url = rule.get("denovo_url")
    if not denovo_url:
        logger.warning(
            f"No de-novo database is configured for db_version '{rule.get('db_version', 'unknown')}'; "
            "skipping de-novo database download."
        )
        return None

    denovo_db_version = rule.get("denovo_db_version")
    if not denovo_db_version:
        raise RuntimeError(f"DB rule defines 'denovo_url' but is missing 'denovo_db_version': {rule}")

    denovo_tar_name = rule.get("denovo_tar_name") or _tar_name_from_url(denovo_url)
    denovo_tar_path = dest_dir / denovo_tar_name

    denovo_extract_dir_name = rule.get("denovo_extract_dir") or _extract_dir_name_from_tar_name(denovo_tar_name)
    denovo_extract_dir = dest_dir / denovo_extract_dir_name

    logger.info(f"Using de-novo database version: {denovo_db_version}")
    logger.info(f"De-novo database source URL: {denovo_url}")

    if denovo_extract_dir.exists() and not force:
        logger.error(f"De-novo database already exists at {denovo_extract_dir}. Use --force to overwrite.")
        raise RuntimeError(f"De-novo database already exists at {denovo_extract_dir}. Use --force to overwrite.")

    _download_file(denovo_url, denovo_tar_path, force=force)
    _extract_tar_gz(denovo_tar_path, denovo_extract_dir, force=force)

    if denovo_tar_path.exists():
        denovo_tar_path.unlink()

    logger.info(f"De-novo database download and extraction complete. Location: {denovo_extract_dir}")
    return denovo_extract_dir

def download_db(
        dest: str,
        checkamg_version: str,
        force: bool = False,
        db_version: str | None = None,
        download_annotate: bool = True,
        download_denovo: bool = True,
    ) -> None:
    if db_version is not None:
        rule = get_db_rule_for_db_version(db_version)
    else:
        rule = get_db_rule_for_checkamg(checkamg_version)

    dest_dir = Path(dest).resolve()
    dest_dir.mkdir(parents=True, exist_ok=True)

    if download_annotate:
        download_annotate_db(rule, dest_dir, checkamg_version, force=force)
    else:
        logger.warning(
            "Skipping annotate database download (--no-download-annotate-db). "
            "The 'checkamg annotate' module will not be usable without it."
        )

    if download_denovo:
        download_denovo_db(rule, dest_dir, force=force)
    else:
        logger.warning(
            "Skipping de-novo database download (--no-download-denovo-db). "
            "The 'checkamg denovo' module will not be usable without it."
        )

    logger.info(f"Database download and extraction complete. Location: {dest_dir}")

def remove_human_readable_files(dest: str) -> None:
    dest_dir = Path(dest).resolve()
    if not dest_dir.exists():
        logger.info(f"Destination does not exist; nothing to remove: {dest_dir}")
        return

    hmm_files = list(dest_dir.rglob("*.hmm"))
    if not hmm_files:
        logger.info(f"No *.hmm files found under: {dest_dir}")
        return

    logger.info(f"Removing {len(hmm_files)} *.hmm files under: {dest_dir}")
    removed = 0
    for fp in hmm_files:
        try:
            fp.unlink()
            removed += 1
        except FileNotFoundError:
            pass
    logger.info(f"Removed {removed} *.hmm files.")

def main(
        command: str,
        dest: str,
        checkamg_version: str,
        force: bool = False,
        db_version: str | None = None,
        rm_hmm: bool = False,
        download_annotate: bool = True,
        download_denovo: bool = True,
    ) -> None:

    logger.info(f"CheckAMG version {checkamg_version}")
    logger.info(f"Starting CheckAMG download...")
    logger.info("Command executed: %s", command)

    download_db(
        dest=dest,
        checkamg_version=checkamg_version,
        force=force,
        db_version=db_version,
        download_annotate=download_annotate,
        download_denovo=download_denovo,
    )

    if rm_hmm:
        remove_human_readable_files(dest=dest)
