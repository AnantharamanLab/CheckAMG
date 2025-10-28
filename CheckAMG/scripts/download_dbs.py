import os
import logging
import requests
import shutil
import gzip
import tarfile
from pathlib import Path
from urllib.request import urlretrieve
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict
from pyhmmer import easel, plan7, hmmer

log_level = logging.INFO
logging.basicConfig(
    level=log_level,
    format='%(asctime)s | %(levelname)s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger()

KEGG_FTP = 'ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz'
KEGG_HTTPS = 'https://www.genome.jp/ftp/db/kofam/profiles.tar.gz'
PFAM_FTP = 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz'
PFAM_HTTPS = 'https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz'
FOAM_FTP = None
FOAM_HTTPS = 'https://osf.io/download/bdpv5'
PHROGS_FTP = None
PHROGS_HTTPS = 'https://phrogs.lmge.uca.fr/downloads_from_website/MSA_phrogs.tar.gz'
DBCAN_FTP = None
DBCAN_HTTPS = 'https://bcb.unl.edu/dbCAN2/download/Databases/V14/dbCAN-HMMdb-V14.txt'
METABOLIC_FTP = None
METABOLIC_HTTPS = 'https://github.com/AnantharamanLab/CheckAMG/raw/refs/heads/main/custom_dbs/METABOLIC_custom.hmm.gz'

def try_download(label, dest, ftp_url, https_url):
    if ftp_url:
        try:
            logger.info(f"Trying FTP download for {label}...")
            urlretrieve(ftp_url, dest)
            logger.info(f"{label} downloaded via FTP.")
        except Exception as ftp_err:
            logger.warning(f"FTP failed for {label}, falling back to HTTPS: {ftp_err}")
            r = requests.get(https_url, stream=True)
            with open(dest, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
            logger.info(f"{label} downloaded via HTTPS.")
    else:
        try:
            logger.info(f"Trying HTTPS download for {label}...")
            r = requests.get(https_url, stream=True)
            with open(dest, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
            logger.info(f"{label} downloaded via HTTPS.")
        except Exception as e:
            raise RuntimeError(f"Failed to download {label}: {e}")

def hmm_db_complete(dest_path):
    prefix = str(dest_path).replace('.hmm', '')
    required = [f"{prefix}.h3m", f"{prefix}.h3i", f"{prefix}.h3f", f"{prefix}.h3p"]
    return all(Path(f).exists() for f in required)

def fix_hmm_names(file):
    # Fix HMM names to be unique by appending a count suffix
    # Even if the order of the HMMs changes in the source file,
    # this shouldn't affect mapping to descritions used by CheckAMG,
    # since the 'ACC" field is used for matching, and those are
    # unique in FOAM.
    logger.info(f"Making HMM names unique in {file}")
    unique_counts = defaultdict(int)
    output = []
    with open(file) as infile:
        block = []
        for line in infile:
            if line.startswith("//"):
                # At end of block, fix name
                for i, l in enumerate(block):
                    if l.startswith("NAME"):
                        name = l.strip().split()[1]
                        unique_counts[name] += 1
                        new_name = f"{name}_{unique_counts[name]}"
                        block[i] = f"NAME  {new_name}\n"
                        break
                output.extend(block + [line])
                block = []
            else:
                block.append(line)
    with open(file, "w") as out:
        out.writelines(output)
        
def build_hmm_from_fasta(msa_path: Path, output_path: Path):
    alphabet = easel.Alphabet.amino()
    builder = plan7.Builder(alphabet)
    background = plan7.Background(alphabet)
    with easel.MSAFile(str(msa_path), digital=True, alphabet=alphabet) as msa_file:
        msa = msa_file.read()
        msa.name = msa.accession = msa_path.stem.encode()
        profile, _, _ = builder.build_msa(msa, background)
        with open(output_path, 'wb') as f:
            profile.write(f)

def build_all_phrog_hmms(msa_dir: Path, out_path: Path, threads: int = 10):
    msa_subdirs = [d for d in msa_dir.iterdir() if d.is_dir()]
    if not msa_subdirs:
        raise RuntimeError("No subdirectory with MSA files found in extracted PHROG archive.")
    msa_data_dir = msa_subdirs[0]
    tmp_hmm_dir = msa_dir / "phrog_hmms"
    tmp_hmm_dir.mkdir(parents=True, exist_ok=True)
    logger.info("Building HMMs from PHROG MSAs...")
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(build_hmm_from_fasta, msa_file, tmp_hmm_dir / (msa_file.stem + ".hmm")) for msa_file in msa_data_dir.glob("*.fma")]
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                logger.error(f"Failed to build HMM: {e}")
    merge_hmm_files_from_dir(tmp_hmm_dir, out_path)
    logger.info(f"Merged PHROG HMMs into {out_path}")
    hmmpress_file(out_path)

def merge_hmm_files_from_dir(src_dir, output_path):
    logger.info(f"Merging HMM files from {src_dir} into {output_path}")
    with open(output_path, 'wb') as out_f:
        for hmm_file in sorted(Path(src_dir).rglob('*.hmm')):
            if hmm_file.is_file() and hmm_file.stat().st_size > 0:
                with open(hmm_file, 'rb') as in_f:
                    shutil.copyfileobj(in_f, out_f)
    shutil.rmtree(src_dir)

def hmmpress_file(hmm_path):
    logger.info(f"Pressing HMM file {hmm_path}")
    hmms = list(plan7.HMMFile(hmm_path))
    output_prefix = str(hmm_path).replace('.hmm', '')
    for ext in ['.h3m', '.h3i', '.h3f', '.h3p']:
        p = Path(f"{output_prefix}{ext}")
        if p.exists():
            p.unlink()
    hmmer.hmmpress(hmms, output_prefix)
    logger.info(f"Pressed HMM database written to {output_prefix}.h3*")

def download_database(label, dest_path, ftp_url, https_url, force=False, decompress=False, untar=False, merge=False, threads=10):
    if hmm_db_complete(dest_path) and not force:
        logger.info(f"{label} already downloaded.")
        return
    tmp_path = dest_path.with_suffix('.tmp')
    try_download(label, tmp_path, ftp_url, https_url)
    if untar:
        extract_dir = dest_path.parent / f"{label}_extracted"
        if extract_dir.exists():
            shutil.rmtree(extract_dir)
        extract_dir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(tmp_path, 'r:gz') as tar:
            tar.extractall(path=extract_dir)
        tmp_path.unlink()
        logger.info(f"{label} downloaded and extracted to {extract_dir}")
        if label == "PHROGs":
            build_all_phrog_hmms(extract_dir, dest_path, threads=threads)
            shutil.rmtree(extract_dir)
        elif merge:
            merge_hmm_files_from_dir(extract_dir, dest_path)
            hmmpress_file(dest_path)
    elif decompress:
        with gzip.open(tmp_path, 'rb') as f_in, open(dest_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        tmp_path.unlink()
        if label == "FOAM":
            fix_hmm_names(dest_path)
        hmmpress_file(dest_path)
    else:
        tmp_path.rename(dest_path)
        hmmpress_file(dest_path)

def remove_human_readable_files(dest):
    to_remove = []
    
    for file in Path(dest).rglob('*.hmm'):
        if not file.name.endswith('.h3i') and not file.name.endswith('.h3m') and not file.name.endswith('.h3p') and not file.name.endswith('.h3f'):
            to_remove.append(file)
            
    if not to_remove:
        logger.info("No human-readable HMM files found to remove.")
        return
    
    logger.info(f"Removing human-readable HMM files from {dest}")
    for file in to_remove:
        os.remove(file)
    logger.info("Human-readable HMM files removed.")

def open_text_auto_gz(path):
    with open(path, 'rb') as probe:
        head = probe.read(2)
    if head == b'\x1f\x8b':
        return gzip.open(path, 'rt')
    return open(path, 'r')

def get_thresholds_from_foam(foam_hmm_path, dest_path):
    logger.info(f"Extracting FOAM thresholds from {foam_hmm_path}")
    acc = None
    cutoff_map = {}
    with open(foam_hmm_path, 'r') as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('ACC'):
                parts = line.split()
                if len(parts) >= 2:
                    acc = parts[1].strip(';')
            elif line.startswith('TC') and acc:
                parts = line.replace(';', '').split()
                full = float(parts[1]) if len(parts) > 1 else None
                dom = float(parts[2]) if len(parts) > 2 else None
                cutoff_map[acc] = (full, dom)
                acc = None
    with open(dest_path, 'w') as out:
        out.write('id\tcutoff_full\tcutoff_domain\n')
        for acc, tup in cutoff_map.items():
            full = '' if tup[0] is None else f"{tup[0]}"
            dom = '' if tup[1] is None else f"{tup[1]}"
            out.write(f"{acc}\t{full}\t{dom}\n")
    logger.info(f"FOAM thresholds written to {dest_path}")

def get_thresholds_from_kegg(ftp_url, https_url, dest_path):
    tmp_path = Path(dest_path).with_suffix('.tmp')
    try_download("KEGG thresholds", tmp_path, ftp_url, https_url)
    with open_text_auto_gz(tmp_path) as f_in, open(dest_path, 'w') as f_out:
        f_out.write('id\tthreshold\n')
        for idx, raw in enumerate(f_in):
            line = raw.rstrip('\n')
            if idx == 0:
                cols = [c.strip() for c in line.split('\t')]
                assert len(cols) >= 2, "ko_list header must have at least two columns"
                assert cols[0].lower() == 'knum', f"Expected first column 'knum', got '{cols[0]}'"
                assert cols[1].lower() == 'threshold', f"Expected second column 'threshold', got '{cols[1]}'"
                continue
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if parts and parts[0].startswith('K'):
                knum = parts[0].strip()
                thr = parts[1].strip() if len(parts) > 1 else ''
                f_out.write(f"{knum}\t{thr}\n")
    try:
        tmp_path.unlink()
    except Exception:
        pass
    logger.info(f"KEGG thresholds written to {dest_path}")
        
def download_all(dest=None, force=False, threads=10):
    os.makedirs(dest, exist_ok=True)
    
    logger.info("Starting download of all databases.")
    dbs = [
        ("KEGG", 'KEGG.hmm', KEGG_FTP, KEGG_HTTPS, True, True, True),
        ("Pfam", 'Pfam-A.hmm', PFAM_FTP, PFAM_HTTPS, True, False, False),
        ("FOAM", 'FOAM.hmm', FOAM_FTP, FOAM_HTTPS, True, False, False),
        ("PHROGs", 'PHROGs.hmm', PHROGS_FTP, PHROGS_HTTPS, False, True, True),
        ("dbCAN", 'dbCAN_HMMdb_v14.hmm', DBCAN_FTP, DBCAN_HTTPS, False, False, False),
        ("METABOLIC", 'METABOLIC_custom.hmm', METABOLIC_FTP, METABOLIC_HTTPS, True, False, False)
    ]
    exceptions = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(download_database, name, Path(dest)/fname, ftp, https, force, decomp, untar, merge, threads): name
            for name, fname, ftp, https, decomp, untar, merge in dbs
        }
        for future in as_completed(futures):
            name = futures[future]
            try:
                future.result()
            except Exception as e:
                logger.error(f"Error downloading {name}: {e}")
                exceptions.append(name)
    if exceptions:
        raise Exception(f"Download failed for: {', '.join(exceptions)}")
    else:
        logger.info("All databases downloaded successfully.")

    # Build thresholds for FOAM and KEGG
    try:
        foam_hmm = Path(dest) / 'FOAM.hmm'
        foam_thr = Path(dest) / 'FOAM_cutoffs.tsv'
        if foam_hmm.exists():
            get_thresholds_from_foam(foam_hmm, foam_thr)
        else:
            logger.warning("FOAM.hmm not found; skipping FOAM thresholds.")
    except Exception as e:
        logger.error(f"FOAM thresholds failed: {e}")
        exceptions.append("FOAM_thresholds")

    try:
        KO_LIST_FTP = KEGG_FTP.replace('profiles.tar.gz', 'ko_list.gz')
        KO_LIST_HTTPS = KEGG_HTTPS.replace('profiles.tar.gz', 'ko_list.gz')
        kegg_thr = Path(dest) / 'KEGG_cutoffs.tsv'
        get_thresholds_from_kegg(KO_LIST_FTP, KO_LIST_HTTPS, kegg_thr)
    except Exception as e:
        logger.error(f"KEGG thresholds failed: {e}")
        exceptions.append("KEGG_thresholds")

    if exceptions:
        raise Exception(f"Completed with errors in: {', '.join(exceptions)}")
    else:
        logger.info("All database thresholds prepared successfully.")