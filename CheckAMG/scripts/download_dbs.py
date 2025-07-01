import os
import requests
from pathlib import Path

def download_pfam(dest):
    url = 'https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam34.0/Pfam-A.hmm.gz'
    dest_path = Path(dest) / 'Pfam-A.hmm.gz'
    if not dest_path.exists():
        r = requests.get(url, stream=True)
        with open(dest_path, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    else:
        print('Pfam already downloaded.')

if __name__ == '__main__':
    download_pfam('databases/')
