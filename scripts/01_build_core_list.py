
from pathlib import Path
import pandas as pd

# Project paths
BASE_DIR = Path(__file__).resolve().parent.parent
SUMMARY_FILE = (
    BASE_DIR
    / "summaryAnnotations"
    / "summary_annotations.tsv"
)

OUTPUT_DIR = BASE_DIR / "data_clean"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Load the summary annotations
df = pd.read_csv(
    SUMMARY_FILE,
    sep="\t",
    dtype=str,
    low_memory=False
)

initial_rows = len(df)
print(f"Initial rows: {initial_rows}")
print("Initial unique rsIDs:",
      df["Variant/Haplotypes"].str.fullmatch(
          r"rs\d+", na=False
      ).sum())

# Filter 1: retain valid single-variant rsIDs
df = df[
    df["Variant/Haplotypes"].str.fullmatch(
        r"rs\d+", na=False
    )
].copy()
after_rs_filter = len(df)

# Filter 2: require gene, drug, and phenotype category
required = ["Gene", "Drug(s)", "Phenotype Category"]
df = df.dropna(subset=required)

for column in required:
    df = df[
        df[column].str.strip().ne("")
    ]

after_required_filter = len(df)

# Rename the variant field for easier processing
df = df.rename(
    columns={"Variant/Haplotypes": "rsID"}
)

# Preserve the source evidence level as a string.
# It is not a ClinVar classification or review-star score.
df = df.rename(
    columns={
        "Level of Evidence": "pharmgkb_evidence_level"
    }
)

# Deduplicate at rsID + gene + drug level
df = df.drop_duplicates(
    subset=["rsID", "Gene", "Drug(s)"]
).copy()
after_dedup = len(df)

# Export the core table
output_file = OUTPUT_DIR / "core_variants.csv"
df.to_csv(output_file, index=False)

# Export unique rsIDs for downstream gnomAD queries
rsids = sorted(df["rsID"].dropna().unique())

with open(
    OUTPUT_DIR / "rsid_list.txt",
    "w",
    encoding="utf-8"
) as file:
    file.write("\n".join(rsids))

# Record the filter funnel
print("\nFILTER FUNNEL")
print(f"Initial rows:                 {initial_rows}")
print(f"After valid rsID filter:      {after_rs_filter}")
print(f"After required fields:        {after_required_filter}")
print(f"After deduplication:          {after_dedup}")
print(f"Unique rsIDs retained:        {len(rsids)}")
print(f"\nSaved: {output_file}")
print(f"Saved: {OUTPUT_DIR / 'rsid_list.txt'}")
