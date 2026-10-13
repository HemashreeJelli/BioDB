
from pathlib import Path
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent
CORE_FILE = BASE_DIR / "data_clean" / "core_variants.csv"

TARGETS = {"rs10012", "rs10123866", "rs10124893"}

df = pd.read_csv(CORE_FILE, dtype=str)

matches = df[df["rsID"].isin(TARGETS)]

columns = [
    col for col in [
        "rsID",
        "Gene",
        "Drug(s)",
        "Phenotype Category",
        "pharmgkb_evidence_level",
        "URL",
    ]
    if col in matches.columns
]

print(matches[columns].to_string(index=False))
print("\nNumber of matching ClinPGx records:", len(matches))
