from pathlib import Path
import pandas as pd
import numpy as np

project_root = Path(__file__).resolve().parent.parent
raw_dir = project_root / "data" / "raw"
processed_dir = project_root / "data" / "processed"
processed_dir.mkdir(parents=True, exist_ok=True)

matrix_file = raw_dir / "GSE197143_series_matrix.txt"
norm_counts_file = raw_dir / "GSE197143_norm_counts.txt"

meta_rows = []
with open(matrix_file, "r", encoding="utf-8", errors="ignore") as f:
    for line in f:
        if line.startswith("!Sample_title"):
            titles = [x.strip('"\n\r ') for x in line.split("\t")[1:]]
        elif line.startswith("!Sample_geo_accession"):
            accessions = [x.strip('"\n\r ') for x in line.split("\t")[1:]]
        elif line.startswith("!Sample_source_name_ch1"):
            tissues = [x.strip('"\n\r ') for x in line.split("\t")[1:]]
        elif line.startswith("!Sample_characteristics_ch1"):
            meta_rows.append([x.strip('"\n\r ') for x in line.split("\t")[1:]])
        elif line.startswith("!series_matrix_table_begin"):
            break

metadata = pd.DataFrame({
    "sample_title": titles,
    "geo_accession": accessions,
    "tissue": tissues,
    "infection_status": [x.split(":")[-1].strip() for x in meta_rows[0]],
    "treatment": [x.split(":")[-1].strip() for x in meta_rows[2]]
})

metadata["group"] = metadata["infection_status"].apply(lambda x: "Infected" if "influenza" in x.lower() else "Mock")
metadata["condition"] = metadata["treatment"]

metadata.to_csv(processed_dir / "sample_metadata.csv", index=False)

expr_df = pd.read_csv(norm_counts_file, sep=r"\s+")
expr_df.to_csv(processed_dir / "cleaned_expression_matrix.csv", index=False)
