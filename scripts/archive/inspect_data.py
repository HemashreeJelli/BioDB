
from pathlib import Path
import pandas as pd

folders = [
    Path("summaryAnnotations"),
    Path("variantAnnotations"),
]

for folder in folders:
    for file in folder.glob("*.tsv"):
        print("\n" + "=" * 70)
        print(f"FILE: {file}")

        df = pd.read_csv(file, sep="\t", dtype=str, low_memory=False)

        print("Rows:", len(df))
        print("Columns:", df.columns.tolist())
        print("\nFirst 2 rows:")
        print(df.head(2).to_string(index=False))
