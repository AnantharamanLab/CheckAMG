import os
import logging
import requests
from pathlib import Path
import shutil
from pyhmmer import easel, plan7, hmmer
import gzip
import tarfile
from concurrent.futures import ThreadPoolExecutor, as_completed

log_level = logging.INFO
logging.basicConfig(
    level=log_level,
    format='%(asctime)s | %(levelname)s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger()

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
        futures = [
            executor.submit(build_hmm_from_fasta, msa_file, tmp_hmm_dir / (msa_file.stem + ".hmm"))
            for msa_file in msa_data_dir.glob("*.fma")
        ]
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
    hmmer.hmmpress(hmms, output_prefix)
    logger.info(f"Pressed HMM database written to {output_prefix}.h3*")

def download_and_handle(url, dest_path, label, force, decompress=False, untar=False, merge=False, threads=10):
    tmp_path = dest_path.with_suffix('.tmp')
    if untar or decompress:
        extract_dir = dest_path.parent / f'{label}_extracted'
        should_download = force or (untar and not extract_dir.exists()) or (decompress and not dest_path.exists())
        if should_download:
            logger.info(f'Downloading {label}...')
            r = requests.get(url, stream=True)
            with open(tmp_path, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
            if untar:
                if extract_dir.exists():
                    shutil.rmtree(extract_dir)
                extract_dir.mkdir(parents=True, exist_ok=True)
                with tarfile.open(tmp_path, 'r:gz') as tar:
                    tar.extractall(path=extract_dir)
                tmp_path.unlink()
                logger.info(f'{label} downloaded and extracted to {extract_dir}')
                if label == "PHROGs":
                    build_all_phrog_hmms(extract_dir, dest_path, threads=threads)
                    shutil.rmtree(extract_dir)
                elif merge:
                    merge_hmm_files_from_dir(extract_dir, dest_path)
                    hmmpress_file(dest_path)
                    logger.info(f'{label} individual HMMs merged and pressed successfully.')
            elif decompress:
                with gzip.open(tmp_path, 'rb') as f_in, open(dest_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                tmp_path.unlink()
                hmmpress_file(dest_path)
                logger.info(f'{label} downloaded, decompressed, and pressed successfully.')
        else:
            logger.info(f'{label} already downloaded.')
    else:
        if not dest_path.exists() or force:
            logger.info(f'Downloading {label}...')
            r = requests.get(url, stream=True)
            with open(tmp_path, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
            tmp_path.rename(dest_path)
            hmmpress_file(dest_path)
            logger.info(f'{label} downloaded and pressed successfully.')
        else:
            logger.info(f'{label} already downloaded.')

def download_kegg(dest, force, threads=10):
    download_and_handle('https://www.genome.jp/ftp/db/kofam/profiles.tar.gz', Path(dest)/'KEGG.hmm', 'KEGG', force, untar=True, merge=True, threads=threads)

def download_pfam(dest, force, threads=10):
    download_and_handle('https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz', Path(dest)/'Pfam-A.hmm', 'Pfam', force, decompress=True, threads=threads)

def download_foam(dest, force, threads=10):
    download_and_handle('https://osf.io/download/bdpv5', Path(dest)/'FOAM.hmm', 'FOAM', force, decompress=True, threads=threads)

def download_phrog(dest, force, threads=10):
    download_and_handle('https://phrogs.lmge.uca.fr/downloads_from_website/MSA_phrogs.tar.gz', Path(dest)/'PHROGs.hmm', 'PHROGs', force, untar=True, threads=threads)

def download_dbcan(dest, force, threads=10):
    download_and_handle('https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-old@UGA/dbCAN-HMMdb-V13.txt', Path(dest)/'dbCAN_HMMdb_v13.hmm', 'dbCAN', force, threads=threads)

def download_metabolic(dest, force, threads=10):
    download_and_handle('https://github.com/AnantharamanLab/CheckAMG/raw/refs/heads/main/custom_dbs/METABOLIC_custom.hmm.gz', Path(dest)/'METABOLIC_custom.hmm', 'METABOLIC', force, decompress=True, threads=threads)

def download_all(dest=None, force=False, threads=10):
    logger.info('Starting download of all databases. This will take ~1 hour depending on your internet connection (KEGG is large and slow to download).')
    os.makedirs(dest, exist_ok=True)
    funcs = [
        ('KEGG', download_kegg),
        ('Pfam', download_pfam),
        ('FOAM', download_foam),
        ('PHROG', download_phrog),
        ('dbCAN', download_dbcan),
        ('METABOLIC', download_metabolic)
    ]
    exceptions = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(func, dest, force, threads): name for name, func in funcs}
        for future in as_completed(futures):
            name = futures[future]
            try:
                future.result()
            except Exception as e:
                logger.error(f'Error downloading {name}: {e}')
                exceptions.append(name)
    if exceptions:
        logger.error(f'The following databases could not be downloaded: {", ".join(exceptions)}')
        raise Exception(f'Download failed for: {", ".join(exceptions)}')
    else:
        logger.info('All databases downloaded successfully.')
