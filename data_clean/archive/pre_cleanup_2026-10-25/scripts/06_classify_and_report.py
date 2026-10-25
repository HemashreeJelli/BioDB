"""
06_classify_and_report.py
=========================
Step 6 of the pharmacogenomic pipeline: 4-group classification.

Input:  data_clean/merged_table.csv

Outputs:
  data_clean/final_table.csv         -- merged_table + group labels + evidence flags
  results/group_counts.csv           -- group counts (rsID level and row level)
  results/top_genes_by_group.csv     -- top genes in each group
  results/top_drugs_by_group.csv     -- top drugs in each group

Classification rules (locked before inspecting counts):
─────────────────────────────────────────────────────────
Skew dimension (primary threshold = 5×):
  high_skew : skew_assessable=True  AND  skew_flag_5x=True
  low_skew  : skew_assessable=True  AND  skew_flag_5x=False
  not_assessable : skew_assessable=False

Evidence dimension — TWO independent views kept separate:

  View A — PharmGKB evidence level (pharmacogenomic-specific):
    high_evidence : level in {1A, 1B, 2A, 2B}
    low_evidence  : level in {3, 4}

  View B — ClinVar review strength (clinical assertion strength):
    high_evidence : clinvar_max_stars >= 3  (expert panel / practice guideline)
    low_evidence  : clinvar_max_stars <= 2  OR  not in ClinVar
    Note: ClinVar stars measure CLINICAL ASSERTION STRENGTH, not
          pharmacogenomic evidence quality. Keep these views SEPARATE.

Primary group assignment uses PharmGKB (View A).
A parallel ClinVar-based group column is also computed (View B).

Four groups (PharmGKB-based, primary):
  Group 1 — HIGH SKEW / HIGH PGKB EVIDENCE  : well-supported disparity
  Group 2 — HIGH SKEW / LOW  PGKB EVIDENCE  : HEADLINE evidence-gap group
  Group 3 — LOW  SKEW / HIGH PGKB EVIDENCE  : reliable evidence, small disparity
  Group 4 — LOW  SKEW / LOW  PGKB EVIDENCE  : lower priority
  Group X — NOT ASSESSABLE                  : insufficient population data

Unit of analysis: rows = rsID × gene × drug relationships.
Unique rsID counts are reported separately to avoid inflation.
"""

import pandas as pd
import numpy as np
from pathlib import Path

ROOT    = Path(__file__).resolve().parent.parent
CLEAN   = ROOT / "data_clean"
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

# ── Load ────────────────────────────────────────────────────────────────────
df = pd.read_csv(CLEAN / "merged_table.csv")
print(f"Loaded merged_table: {len(df):,} rows, {df['rsID'].nunique():,} unique rsIDs")

# ── Evidence rules (defined before inspecting output) ───────────────────────
HIGH_PGKB = {"1A", "1B", "2A", "2B"}
HIGH_CLINVAR_STARS = 3   # >= 3 stars = expert panel or practice guideline

# ── Classify PharmGKB evidence ──────────────────────────────────────────────
df["pgkb_evidence_tier"] = df["pharmgkb_evidence_level"].apply(
    lambda x: "high" if str(x).strip() in HIGH_PGKB else "low"
)

# ── Classify ClinVar evidence ───────────────────────────────────────────────
def clinvar_tier(row):
    if not row["in_clinvar"]:
        return "not_in_clinvar"
    stars = row["clinvar_max_stars"]
    if pd.isna(stars):
        return "not_in_clinvar"
    return "high" if stars >= HIGH_CLINVAR_STARS else "low"

df["clinvar_evidence_tier"] = df.apply(clinvar_tier, axis=1)

# ── Classify skew ───────────────────────────────────────────────────────────
def skew_tier(row):
    if not row["skew_assessable"]:
        return "not_assessable"
    return "high" if row["skew_flag_5x"] else "low"

df["skew_tier"] = df.apply(skew_tier, axis=1)

# ── Assign primary group (PharmGKB-based) ───────────────────────────────────
GROUP_LABELS = {
    ("high", "high"): "Group 1: High Skew / High Evidence",
    ("high", "low"):  "Group 2: High Skew / Low Evidence",
    ("low",  "high"): "Group 3: Low Skew / High Evidence",
    ("low",  "low"):  "Group 4: Low Skew / Low Evidence",
}

def assign_pgkb_group(row):
    if row["skew_tier"] == "not_assessable":
        return "Group X: Not Assessable"
    key = (row["skew_tier"], row["pgkb_evidence_tier"])
    return GROUP_LABELS[key]

df["group_pgkb"] = df.apply(assign_pgkb_group, axis=1)

# ── Assign parallel ClinVar-based group ─────────────────────────────────────
def assign_clinvar_group(row):
    if row["skew_tier"] == "not_assessable":
        return "Group X: Not Assessable"
    cv_tier = row["clinvar_evidence_tier"]
    if cv_tier == "not_in_clinvar":
        cv_side = "low"
    else:
        cv_side = cv_tier   # "high" or "low"
    key = (row["skew_tier"], cv_side)
    return GROUP_LABELS[key]

df["group_clinvar"] = df.apply(assign_clinvar_group, axis=1)

# ── Save final table ─────────────────────────────────────────────────────────
out_final = CLEAN / "final_table.csv"
df.to_csv(out_final, index=False)
print(f"Saved final_table.csv: {df.shape}")

# ── Group counts ─────────────────────────────────────────────────────────────
print("\n" + "="*60)
print("GROUP COUNTS — PharmGKB-based (primary)")
print("="*60)
print("\nUnit: rsID × gene × drug ROWS")
row_counts_pgkb = df["group_pgkb"].value_counts().sort_index()
for group, count in row_counts_pgkb.items():
    pct = 100 * count / len(df)
    print(f"  {group:<45} {count:>5}  ({pct:.1f}%)")

print("\nUnit: unique rsIDs")
rsid_counts_pgkb = df.groupby("group_pgkb")["rsID"].nunique().sort_index()
for group, count in rsid_counts_pgkb.items():
    pct = 100 * count / df["rsID"].nunique()
    print(f"  {group:<45} {count:>5}  ({pct:.1f}%)")

print("\n" + "="*60)
print("GROUP COUNTS — ClinVar-based (secondary / parallel view)")
print("="*60)
print("\nUnit: rsID × gene × drug ROWS")
row_counts_cv = df["group_clinvar"].value_counts().sort_index()
for group, count in row_counts_cv.items():
    pct = 100 * count / len(df)
    print(f"  {group:<45} {count:>5}  ({pct:.1f}%)")

# ── Save group_counts.csv ───────────────────────────────────────────────────
gc = pd.DataFrame({
    "group_pgkb":           row_counts_pgkb.index,
    "rows_pgkb":            row_counts_pgkb.values,
    "pct_rows_pgkb":        (100 * row_counts_pgkb / len(df)).round(1).values,
    "unique_rsids_pgkb":    rsid_counts_pgkb.reindex(row_counts_pgkb.index).values,
})
cv_row = row_counts_cv.reindex(gc["group_pgkb"]).values
gc["rows_clinvar_view"] = cv_row
gc.to_csv(RESULTS / "group_counts.csv", index=False)
print(f"\nSaved: results/group_counts.csv")

# ── Top genes & drugs in each group (PharmGKB) ──────────────────────────────
print("\n" + "="*60)
print("TOP GENES BY GROUP (PharmGKB, top 10 per group)")
print("="*60)
gene_rows = []
for group in sorted(df["group_pgkb"].unique()):
    sub = df[df["group_pgkb"] == group]
    top = sub.groupby("Gene")["rsID"].nunique().sort_values(ascending=False).head(10)
    for gene, n in top.items():
        gene_rows.append({"group": group, "gene": gene, "unique_rsids": n})
    print(f"\n  {group}")
    print("  " + ", ".join(f"{g}({n})" for g, n in top.items()))

gene_df = pd.DataFrame(gene_rows)
gene_df.to_csv(RESULTS / "top_genes_by_group.csv", index=False)
print(f"\nSaved: results/top_genes_by_group.csv")

print("\n" + "="*60)
print("TOP DRUGS BY GROUP (PharmGKB, top 10 per group)")
print("="*60)
drug_rows = []
for group in sorted(df["group_pgkb"].unique()):
    sub = df[df["group_pgkb"] == group]
    top = sub["Drug(s)"].value_counts().head(10)
    for drug, n in top.items():
        drug_rows.append({"group": group, "drug": drug, "annotation_rows": n})
    print(f"\n  {group}")
    for drug, n in top.items():
        print(f"    {drug} ({n})")

drug_df = pd.DataFrame(drug_rows)
drug_df.to_csv(RESULTS / "top_drugs_by_group.csv", index=False)
print(f"\nSaved: results/top_drugs_by_group.csv")

# ── Headline Group 2 deep look ───────────────────────────────────────────────
g2 = df[df["group_pgkb"] == "Group 2: High Skew / Low Evidence"]
print("\n" + "="*60)
print("HEADLINE GROUP 2 — High Skew / Low PharmGKB Evidence")
print("="*60)
print(f"  Rows:        {len(g2):,}")
print(f"  Unique rsIDs:{g2['rsID'].nunique():,}")
print(f"  in_clinvar:  {g2['in_clinvar'].sum():,} ({100*g2['in_clinvar'].mean():.1f}%)")
print(f"  zero_min_af: {g2['zero_min_af_flag'].sum():,}")
print(f"\n  ClinVar star distribution within Group 2:")
print(g2["clinvar_evidence_tier"].value_counts(dropna=False).to_string())
print(f"\n  Phenotype categories:")
print(g2["Phenotype Category"].value_counts().to_string())
print(f"\n  Top 10 genes:")
print(g2.groupby("Gene")["rsID"].nunique().sort_values(ascending=False).head(10).to_string())
print(f"\n  Top 10 highest disparity_ratio variants:")
top_g2 = g2.nlargest(10, "disparity_ratio")[
    ["rsID", "Gene", "Drug(s)", "max_pop", "max_AF",
     "min_pop", "min_AF", "disparity_ratio", "absolute_AF_diff",
     "pharmgkb_evidence_level", "clinvar_evidence_tier"]
]
print(top_g2.to_string(index=False))

# ── Sensitivity check: 3x and 10x thresholds ────────────────────────────────
print("\n" + "="*60)
print("SENSITIVITY CHECK — Group counts at 3× and 10× thresholds")
print("="*60)
for thresh, col in [("3x", "skew_flag_3x"), ("10x", "skew_flag_10x")]:
    def assign_sensitivity(row, flag_col=col):
        if not row["skew_assessable"]:
            return "Group X: Not Assessable"
        skew = "high" if row[flag_col] else "low"
        key = (skew, row["pgkb_evidence_tier"])
        return GROUP_LABELS[key]
    tmp = df.apply(assign_sensitivity, axis=1).value_counts().sort_index()
    print(f"\n  Threshold: {thresh}")
    for g, n in tmp.items():
        print(f"    {g:<45} {n:>5}")

print("\nDone.")
