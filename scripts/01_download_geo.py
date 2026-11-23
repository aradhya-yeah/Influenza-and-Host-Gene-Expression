from pathlib import Path
import urllib.request
import gzip
import shutil

project_root = Path(__file__).resolve().parent.parent
raw_dir = project_root / "data" / "raw"
raw_dir.mkdir(parents=True, exist_ok=True)

matrix_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE197nnn/GSE197143/matrix/GSE197143_series_matrix.txt.gz"
suppl_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE197nnn/GSE197143/suppl/GSE197143_btchCORR_norm_Brunotte_Hulu_FluB_2022_120122.txt.gz"

matrix_gz = raw_dir / "GSE197143_series_matrix.txt.gz"
matrix_txt = raw_dir / "GSE197143_series_matrix.txt"
suppl_gz = raw_dir / "GSE197143_norm_counts.txt.gz"
suppl_txt = raw_dir / "GSE197143_norm_counts.txt"

headers = {"User-Agent": "InfluenzaTranscriptomics-UndergradResearch/1.0"}

if not matrix_gz.exists():
    req = urllib.request.Request(matrix_url, headers=headers)
    with urllib.request.urlopen(req, timeout=60) as response, open(matrix_gz, "wb") as out_file:
        shutil.copyfileobj(response, out_file)

if not matrix_txt.exists():
    with gzip.open(matrix_gz, "rb") as f_in, open(matrix_txt, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

if not suppl_gz.exists():
    req = urllib.request.Request(suppl_url, headers=headers)
    with urllib.request.urlopen(req, timeout=120) as response, open(suppl_gz, "wb") as out_file:
        shutil.copyfileobj(response, out_file)

if not suppl_txt.exists():
    with gzip.open(suppl_gz, "rb") as f_in, open(suppl_txt, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
