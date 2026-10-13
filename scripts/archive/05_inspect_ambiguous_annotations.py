
from pathlib import Path
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent

TARGETS = {"rs10012", "rs10123866", "rs10124893"}

FILES = [
    BASE_DIR / "summaryAnnotations" / "summary_annotations.tsv",
    BASE_DIR / "variantAnnotations" / "var_drug_ann.tsv",
    BASE_DIR / "variantAnnotations" / "var_pheno_ann.tsv",
]

for file_path in FILES:
    if not file_path.exists():
        print(f"\nFile not found: {file_path}")
        continue

    print(f"\n{'=' * 70}")
    print(f"FILE: {file_path.name}")

    df = pd.read_csv(file_path, sep="\t", dtype=str)

    if "Variant/Haplotypes" not in df.columns:
        print("No Variant/Haplotypes column.")
        continue

    mask = df["Variant/Haplotypes"].fillna("").apply(
        lambda value: any(
            target in value.split(",") or target == value.strip()
            for target in TARGETS
        )
    )

    matches = df.loc[mask]
    print("Matching rows:", len(matches))

    if not matches.empty:
        print(matches.to_string(index=False, max_cols=None))
