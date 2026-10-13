
from pathlib import Path
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent

TARGET_IDS = {
    "1450928473",  # rs10012
    "1449168475",  # rs10123866
    "1449168489",  # rs10124893
}

FILE = (
    BASE_DIR
    / "summaryAnnotations"
    / "summary_ann_evidence.tsv"
)

df = pd.read_csv(FILE, sep="\t", dtype=str)

matches = df[df["Summary Annotation ID"].isin(TARGET_IDS)]

print(matches.to_string(index=False, max_cols=None))
print("\nMatching evidence records:", len(matches))
