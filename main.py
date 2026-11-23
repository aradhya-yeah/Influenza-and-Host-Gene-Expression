import subprocess
import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent
scripts_dir = project_root / "scripts"

subprocess.run([sys.executable, str(scripts_dir / "01_download_geo.py")], check=True)
subprocess.run([sys.executable, str(scripts_dir / "02_clean_and_inspect.py")], check=True)
subprocess.run([sys.executable, str(scripts_dir / "03_differential_expression.py")], check=True)
subprocess.run([sys.executable, str(scripts_dir / "04_generate_figures.py")], check=True)
