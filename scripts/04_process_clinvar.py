"""
04_process_clinvar.py
=====================
Step 3 of the pharmacogenomic pipeline.

Input  : variant_summary.txt.gz  (ClinVar bulk download, project root)
         data_clean/rsid_list.txt (2,604 rsIDs from PharmGKB core list)

Output : data_clean/clinvar_matched.csv
           - One aggregated row per rsID (with summary fields)
         data_clean/clinvar_all_records.csv
           - All individual ClinVar records for each matched rsID (multi-row)
           - Useful for inspecting conflicts / condition-specific assertions

Decisions locked before running:
  - Only GRCh38 records are used.
  - rsIDs with RS# == -1 (no dbSNP link) are excluded.
  - ReviewStatus → star score mapping is explicit and transparent;
    the original ReviewStatus text is always kept.
  - Multiple ClinVar records per rsID are preserved in the _all_records file
    and AGGREGATED (not reduced to just the highest star) in the _matched file.
  - 'drug response' ClinicalSignificance is flagged separately.
  - rsIDs not found in ClinVar are flagged with in_clinvar = False.
"""

import gzip
import pandas as pd
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parent.parent
CLINVAR_GZ   = ROOT / "variant_summary.txt.gz"
RSID_LIST    = ROOT / "data_clean" / "rsid_list.txt"
OUT_MATCHED  = ROOT / "data_clean" / "clinvar_matched.csv"
OUT_ALL      = ROOT / "data_clean" / "clinvar_all_records.csv"

# ── Star-rating mapping (based on actual values inspected 2026-10-12) ──────
# Retain all values seen; anything missing gets starred as 0 with a warning.
STAR_MAP = {
    "practice guideline":                                    4,
    "reviewed by expert panel":                              3,
    "criteria provided, multiple submitters, no conflicts":  2,
    "criteria provided, single submitter":                   1,
    "criteria provided, conflicting classifications":        1,
    "no assertion criteria provided":                        0,
    "no classification provided":                            0,
    "no classification for the single variant":              0,
    "no classifications from unflagged records":             0,
    "-":                                                     0,
}

# ── Load rsID target list ───────────────────────────────────────────────────
target_rsids = set(pd.read_csv(RSID_LIST, header=None)[0].astype(str).str.strip())
print(f"Target rsIDs loaded: {len(target_rsids)}")

# ── Read ClinVar (only needed columns) ─────────────────────────────────────
USECOLS = [
    "RS# (dbSNP)", "GeneSymbol", "Type",
    "ClinicalSignificance", "ClinSigSimple",
    "ReviewStatus", "NumberSubmitters",
    "Assembly", "PhenotypeList",
]

print("Reading ClinVar variant_summary.txt.gz …")
cv = pd.read_csv(CLINVAR_GZ, sep="\t", usecols=USECOLS, low_memory=False)
print(f"  Total rows: {len(cv):,}")

# ── Filter: GRCh38 only ─────────────────────────────────────────────────────
cv = cv[cv["Assembly"] == "GRCh38"].copy()
print(f"  GRCh38 rows: {len(cv):,}")

# ── Filter: valid RS number ─────────────────────────────────────────────────
cv = cv[cv["RS# (dbSNP)"].notna() & (cv["RS# (dbSNP)"] != -1)].copy()
cv["rsID"] = "rs" + cv["RS# (dbSNP)"].astype(int).astype(str)
print(f"  Rows with valid rsID: {len(cv):,}")

# ── Filter: restrict to target rsIDs only ──────────────────────────────────
cv = cv[cv["rsID"].isin(target_rsids)].copy()
print(f"  Rows matching target rsIDs: {len(cv):,}")
print(f"  Unique target rsIDs found in ClinVar: {cv['rsID'].nunique():,}")

# ── Check for unseen ReviewStatus values (safety net) ──────────────────────
unseen = set(cv["ReviewStatus"].unique()) - set(STAR_MAP.keys())
if unseen:
    print(f"\n  WARNING: Unseen ReviewStatus values (will be scored 0): {unseen}")

# ── Derive star score ───────────────────────────────────────────────────────
cv["clinvar_stars"] = cv["ReviewStatus"].map(STAR_MAP).fillna(0).astype(int)

# ── Flag drug response records ─────────────────────────────────────────────
cv["is_drug_response"] = cv["ClinicalSignificance"].str.lower().str.contains(
    "drug response", na=False
)

# ── Save all individual records (pre-aggregation) ──────────────────────────
all_records = cv[[
    "rsID", "GeneSymbol", "Type", "ClinicalSignificance",
    "ReviewStatus", "clinvar_stars", "NumberSubmitters",
    "PhenotypeList", "is_drug_response",
]].rename(columns={"GeneSymbol": "clinvar_gene", "Type": "variant_type"})
all_records.to_csv(OUT_ALL, index=False)
print(f"\nSaved all individual records: {OUT_ALL}")

# ── Aggregate: one row per rsID ─────────────────────────────────────────────
def aggregate_rsid(grp):
    """Summarise all ClinVar records for one rsID."""
    classifications = grp["ClinicalSignificance"].dropna().unique().tolist()
    review_statuses = grp["ReviewStatus"].dropna().unique().tolist()
    max_stars       = int(grp["clinvar_stars"].max())
    n_records       = len(grp)
    n_submitters    = int(grp["NumberSubmitters"].max())  # highest submitter count seen
    has_drug_resp   = grp["is_drug_response"].any()
    has_pathogenic  = grp["ClinicalSignificance"].str.lower().str.contains(
        "pathogenic", na=False).any()
    has_conflict    = grp["ClinicalSignificance"].str.lower().str.contains(
        "conflicting", na=False).any()
    has_uncertain   = grp["ClinicalSignificance"].str.lower().str.contains(
        "uncertain", na=False).any()

    # Highest-star ReviewStatus text (for reference)
    best_row = grp.loc[grp["clinvar_stars"].idxmax()]

    return pd.Series({
        "clinvar_n_records":           n_records,
        "clinvar_max_stars":           max_stars,
        "clinvar_best_review_status":  best_row["ReviewStatus"],
        "clinvar_classification_list": "; ".join(sorted(set(classifications))),
        "clinvar_review_status_list":  "; ".join(sorted(set(review_statuses))),
        "clinvar_max_submitters":      n_submitters,
        "clinvar_has_drug_response":   has_drug_resp,
        "clinvar_has_pathogenic":      has_pathogenic,
        "clinvar_has_conflicting":     has_conflict,
        "clinvar_has_uncertain":       has_uncertain,
        # Primary classification = that of the best-star row
        "clinvar_primary_classification": best_row["ClinicalSignificance"],
    })

print("\nAggregating by rsID …")
agg = cv.groupby("rsID").apply(aggregate_rsid).reset_index()
print(f"  Aggregated rows: {len(agg):,}")

# ── Build final matched table (all target rsIDs, left-join) ────────────────
rsid_df = pd.DataFrame({"rsID": sorted(target_rsids)})
matched = rsid_df.merge(agg, on="rsID", how="left")
matched["in_clinvar"] = matched["clinvar_max_stars"].notna()

# Fill summary flags for rsIDs absent from ClinVar
bool_cols = [
    "clinvar_has_drug_response", "clinvar_has_pathogenic",
    "clinvar_has_conflicting",   "clinvar_has_uncertain",
]
for col in bool_cols:
    matched[col] = matched[col].fillna(False)
matched["clinvar_n_records"] = matched["clinvar_n_records"].fillna(0).astype(int)

# ── Coverage report ─────────────────────────────────────────────────────────
n_total   = len(matched)
n_found   = matched["in_clinvar"].sum()
n_missing = n_total - n_found
print(f"\n=== ClinVar Coverage ===")
print(f"  Target rsIDs:          {n_total:,}")
print(f"  Found in ClinVar:      {n_found:,} ({100*n_found/n_total:.1f}%)")
print(f"  Not in ClinVar:        {n_missing:,} ({100*n_missing/n_total:.1f}%)")
print(f"\n  Star distribution (best star per rsID):")
print(matched[matched["in_clinvar"]]["clinvar_max_stars"].value_counts().sort_index())
print(f"\n  Has drug response flag: {matched['clinvar_has_drug_response'].sum():,}")
print(f"  Has conflicting:        {matched['clinvar_has_conflicting'].sum():,}")

# ── Save matched table ──────────────────────────────────────────────────────
matched.to_csv(OUT_MATCHED, index=False)
print(f"\nSaved: {OUT_MATCHED}")
print("Done.")
