"""
01_build_core_list.py -- Step 1: PharmGKB core list

Input : summaryAnnotations/summary_annotations.tsv
Output: data_clean/core_variants.csv   (one row per rsID x gene x drug)
        data_clean/rsid_list.txt       (unique rsIDs for the gnomAD and ClinVar steps)

Keeps single-rsID annotations (no star alleles / haplotypes) that have a gene,
drug and phenotype category. Duplicates keep the strongest evidence level.
"""

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

# Deduplicate at rsID + gene + drug level, keeping the strongest
# PharmGKB evidence level when the same relationship appears more than once.
LEVEL_RANK = {"1A": 0, "1B": 1, "2A": 2, "2B": 3, "3": 4, "4": 5}
df["_rank"] = df["pharmgkb_evidence_level"].str.strip().map(LEVEL_RANK).fillna(99)
df = df.sort_values("_rank", kind="stable").drop(columns="_rank")
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
