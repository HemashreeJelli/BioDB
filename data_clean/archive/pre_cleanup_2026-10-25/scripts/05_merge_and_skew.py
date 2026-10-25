"""
05_merge_and_skew.py
====================
Steps 4 & 5 of the pharmacogenomic pipeline.

Inputs:
  data_clean/core_rows_usable_v3.csv       -- 4016 core PharmGKB relationships
  data_clean/annotation_resolution_v3.csv  -- per-annotation allele resolution
  data_clean/site_allele_frequencies.csv   -- gnomAD ALT-allele population counts
  data_clean/clinvar_matched.csv           -- ClinVar aggregated per rsID

Output:
  data_clean/merged_table.csv              -- Full merged table (one row per core relationship)
  data_clean/sanity_checks.txt             -- All validation checks, printed and saved

Design decisions (locked):
  - Only 'usable' core rows (3834) carry forward through the merge.
  - For multi-ALT rsIDs, frequency of the pharmacogenomically annotated allele
    is selected using effect_alleles_plus from annotation_resolution_v3.
    If the match is ambiguous, the highest-AC ALT allele is used and flagged.
  - Populations used for skew: AFR, AMR, EAS, NFE, SAS.
    FIN and MID are computed but excluded from primary max/min comparison
    (too small or inconsistently available).
  - Minimum AN threshold: 1000 per population. Populations below this are
    masked (treated as missing) for that variant's skew computation.
  - Zero AF: NOT treated as zero in ratio calculation.
    If min_AF == 0 the ratio is left as NaN and zero_min_af_flag = True.
  - Minimum max_AF filter (secondary): 0.01 (1%). Variants with max_AF < 1%
    across all primary populations get skew_assessable = False.
  - Ratio thresholds flagged: 3x, 5x, 10x.
  - ClinVar and gnomAD joins are LEFT joins -- every usable PharmGKB row is kept.
"""

import pandas as pd
import numpy as np
from pathlib import Path

ROOT      = Path(__file__).resolve().parent.parent
CLEAN     = ROOT / "data_clean"

# ── Config ──────────────────────────────────────────────────────────────────
PRIMARY_POPS  = ["AFR", "AMR", "EAS", "NFE", "SAS"]
MIN_AN        = 1_000          # mask population if AN < this
MIN_MAX_AF    = 0.01           # secondary filter: skip skew if max_AF < 1%

# ── Load inputs ─────────────────────────────────────────────────────────────
print("Loading inputs…")
core  = pd.read_csv(CLEAN / "core_rows_usable_v3.csv")
res   = pd.read_csv(CLEAN / "annotation_resolution_v3.csv")
freq  = pd.read_csv(CLEAN / "site_allele_frequencies.csv")
clin  = pd.read_csv(CLEAN / "clinvar_matched.csv")

print(f"  core rows:              {len(core):,}  (rsIDs: {core['rsID'].nunique():,})")
print(f"  annotation resolutions: {len(res):,}")
print(f"  freq rows (all):        {len(freq):,}  (rsIDs: {freq['rsID'].nunique():,})")
print(f"  clinvar rows:           {len(clin):,}")

# ── Step 1: restrict core to usable rows ────────────────────────────────────
core_usable = core[core["usability"] == "usable"].copy()
print(f"\nUsable core rows: {len(core_usable):,}")

# ── Step 2: build a per-rsID allele lookup from annotation_resolution_v3 ────
# For each rsID, what is the set of accepted effect alleles (plus strand)?
accepted_res = res[res["resolution_v3"].str.startswith("accepted", na=False)].copy()
accepted_res = accepted_res.dropna(subset=["rsID", "effect_alleles_plus"])

# Aggregate accepted effect alleles per rsID (take the most common one if multiple)
def pick_effect_allele(group):
    """Return the single-letter effect allele most frequently cited per rsID."""
    alleles = group["effect_alleles_plus"].dropna()
    # Flatten multi-letter tokens e.g. 'AG' -> {'A','G'}
    letter_sets = alleles.apply(lambda x: set(str(x).upper().replace(" ", "")))
    flat = [l for s in letter_sets for l in s]
    if not flat:
        return None
    # Most common single letter
    from collections import Counter
    counts = Counter(flat)
    return counts.most_common(1)[0][0]

rsid_effect_allele = (
    accepted_res.groupby("rsID")
    .apply(pick_effect_allele)
    .rename("effect_allele_plus")
    .reset_index()
)
print(f"\nrsIDs with resolved effect allele: {len(rsid_effect_allele):,}")

# ── Step 3: build the frequency lookup (ALT rows only) ──────────────────────
alt_freq = freq[freq["allele_type"] == "ALT"].copy()
print(f"ALT frequency rows: {len(alt_freq):,}")

# Identify multi-ALT rsIDs
alt_per_rsid = alt_freq.groupby("rsID").size()
multi_alt_rsids = set(alt_per_rsid[alt_per_rsid > 1].index)
print(f"rsIDs with >1 ALT row: {len(multi_alt_rsids):,}")


print("\nSelecting best ALT row per rsID (vectorized)…")

# Merge effect allele lookup into all ALT rows once
alt_freq_with_effect = alt_freq.merge(rsid_effect_allele, on="rsID", how="left")

# Single-ALT rsIDs: take directly
single_rows = alt_freq_with_effect[~alt_freq_with_effect["rsID"].isin(multi_alt_rsids)].copy()
single_rows["allele_selection_method"] = "single_alt"

# Multi-ALT rsIDs — Attempt 1: exact match on effect allele letter
multi_rows = alt_freq_with_effect[alt_freq_with_effect["rsID"].isin(multi_alt_rsids)].copy()
effect_match = multi_rows[
    multi_rows["allele"].notna() &
    multi_rows["effect_allele_plus"].notna() &
    (multi_rows["allele"] == multi_rows["effect_allele_plus"])
].copy()
effect_match = (
    effect_match.sort_values("overall_ac", ascending=False)
    .drop_duplicates("rsID", keep="first")
    .copy()
)
effect_match["allele_selection_method"] = "effect_allele_match"

# Multi-ALT rsIDs — Attempt 2: fallback to highest-AC for remaining rsIDs
resolved_rsids = set(effect_match["rsID"])
fallback_rows = multi_rows[~multi_rows["rsID"].isin(resolved_rsids)].copy()
fallback_rows = (
    fallback_rows.sort_values("overall_ac", ascending=False)
    .drop_duplicates("rsID", keep="first")
    .copy()
)
fallback_rows["allele_selection_method"] = "highest_ac_fallback"

best_alt = pd.concat([single_rows, effect_match, fallback_rows], ignore_index=True)
print(f"Best ALT rows selected: {len(best_alt):,}")
print("\nAllele selection method:")
print(best_alt["allele_selection_method"].value_counts())

# ── Step 4: Mask populations below AN threshold ──────────────────────────────
print(f"\nMasking populations with AN < {MIN_AN}…")
for pop in PRIMARY_POPS:
    an_col = f"{pop}_an"
    af_col = f"{pop}_af"
    low_mask = best_alt[an_col] < MIN_AN
    if low_mask.sum() > 0:
        print(f"  {pop}: masking {low_mask.sum()} rows (AN < {MIN_AN})")
        best_alt.loc[low_mask, af_col] = np.nan

# ── Step 5: Calculate skew metrics across primary populations ────────────────
print("\nCalculating skew metrics…")
af_cols = [f"{p}_af" for p in PRIMARY_POPS]
an_cols = [f"{p}_an" for p in PRIMARY_POPS]

def compute_skew(row):
    afs = {pop: row[f"{pop}_af"] for pop in PRIMARY_POPS
           if pd.notna(row[f"{pop}_af"]) and row[f"{pop}_an"] >= MIN_AN}

    # Measure skew on the minor allele. If ALT is the major allele on average,
    # use 1 - ALT AF instead, so the result doesn't depend on which allele
    # the reference genome happens to call ALT.
    skew_allele = "ALT"
    if afs and np.mean(list(afs.values())) > 0.5:
        afs = {p: 1 - v for p, v in afs.items()}
        skew_allele = "REF_or_other"

    n_pops_eligible = len(afs)
    if n_pops_eligible < 2:
        return pd.Series({
            "n_pops_eligible": n_pops_eligible,
            "max_AF": np.nan, "min_AF": np.nan,
            "max_pop": np.nan, "min_pop": np.nan,
            "absolute_AF_diff": np.nan, "disparity_ratio": np.nan,
            "zero_min_af_flag": False,
            "skew_flag_3x": False, "skew_flag_5x": False, "skew_flag_10x": False,
            "skew_assessable": False,
            "skew_not_assessable_reason": f"only {n_pops_eligible} eligible population(s)",
        })

    max_pop = max(afs, key=afs.get)
    min_pop = min(afs, key=afs.get)
    max_af  = afs[max_pop]
    min_af  = afs[min_pop]
    abs_diff = max_af - min_af

    # Secondary filter: max_AF < 1% -> not assessable for skew
    if max_af < MIN_MAX_AF:
        return pd.Series({
            "n_pops_eligible": n_pops_eligible,
            "max_AF": max_af, "min_AF": min_af,
            "max_pop": max_pop, "min_pop": min_pop,
            "absolute_AF_diff": abs_diff, "disparity_ratio": np.nan,
            "zero_min_af_flag": min_af == 0,
            "skew_flag_3x": False, "skew_flag_5x": False, "skew_flag_10x": False,
            "skew_assessable": False,
            "skew_not_assessable_reason": f"max_AF {max_af:.4f} < {MIN_MAX_AF} threshold",
        })

    # Ratio
    if min_af == 0:
        ratio = np.nan
        zero_flag = True
    else:
        ratio = max_af / min_af
        zero_flag = False

    # min_AF == 0 with max_AF >= 1% is an unbounded ratio: it passes every threshold.
    return pd.Series({
        "n_pops_eligible": n_pops_eligible,
        "skew_allele": skew_allele,
        "max_AF": max_af, "min_AF": min_af,
        "max_pop": max_pop, "min_pop": min_pop,
        "absolute_AF_diff": abs_diff,
        "disparity_ratio": ratio,
        "zero_min_af_flag": zero_flag,
        "skew_flag_3x":  zero_flag or (ratio >= 3),
        "skew_flag_5x":  zero_flag or (ratio >= 5),
        "skew_flag_10x": zero_flag or (ratio >= 10),
        "skew_assessable": True,
        "skew_not_assessable_reason": None,
    })

skew_cols = best_alt.apply(compute_skew, axis=1)
best_alt = pd.concat([best_alt, skew_cols], axis=1)

# ── Step 6: Merge everything onto the core usable list ──────────────────────
print("\nMerging…")
n_before = len(core_usable)

# Columns to bring from freq
freq_keep = (
    ["rsID", "variant_id", "allele", "allele_selection_method",
     "overall_ac", "overall_an", "overall_af"]
    + [f"{p}_{x}" for p in PRIMARY_POPS + ["FIN","MID"] for x in ["ac","an","af"]]
    + ["n_pops_eligible", "skew_allele", "max_AF", "min_AF", "max_pop", "min_pop",
       "absolute_AF_diff", "disparity_ratio", "zero_min_af_flag",
       "skew_flag_3x", "skew_flag_5x", "skew_flag_10x",
       "skew_assessable", "skew_not_assessable_reason"]
)
freq_keep = [c for c in freq_keep if c in best_alt.columns]
best_alt_slim = best_alt[freq_keep].copy()

merged = core_usable.merge(best_alt_slim, on="rsID", how="left")
n_after = len(merged)

print(f"  Core rows before merge: {n_before:,}")
print(f"  Core rows after merge:  {n_after:,}")
assert n_before == n_after, "ROW COUNT CHANGED after freq merge! Investigate."

# Add gnomAD coverage flag
merged["in_gnomad"] = merged["variant_id"].notna()

# Merge ClinVar
clinvar_keep = [c for c in clin.columns if c not in ["rsID"] or c == "rsID"]
merged = merged.merge(clin[clinvar_keep], on="rsID", how="left")
n_after2 = len(merged)
print(f"  Core rows after ClinVar merge: {n_after2:,}")
assert n_before == n_after2, "ROW COUNT CHANGED after ClinVar merge! Investigate."

# ── Step 7: Sanity checks ────────────────────────────────────────────────────
lines = []
lines.append("=== MERGE & SKEW SANITY CHECKS ===\n")

lines.append(f"Core usable rows (input):       {n_before:,}")
lines.append(f"Merged rows (output):           {n_after2:,}")
lines.append(f"Unique rsIDs (input):           {core_usable['rsID'].nunique():,}")
lines.append(f"Unique rsIDs (merged):          {merged['rsID'].nunique():,}")
lines.append("")

lines.append("--- gnomAD Coverage ---")
lines.append(f"in_gnomad = True:               {merged['in_gnomad'].sum():,} ({100*merged['in_gnomad'].mean():.1f}%)")
lines.append(f"in_gnomad = False:              {(~merged['in_gnomad']).sum():,} ({100*(~merged['in_gnomad']).mean():.1f}%)")
lines.append("")

lines.append("--- ClinVar Coverage ---")
lines.append(f"in_clinvar = True:              {merged['in_clinvar'].sum():,} ({100*merged['in_clinvar'].mean():.1f}%)")
lines.append(f"in_clinvar = False:             {(~merged['in_clinvar']).sum():,} ({100*(~merged['in_clinvar']).mean():.1f}%)")
lines.append("")

lines.append("--- Allele Selection Method ---")
for method, count in merged["allele_selection_method"].value_counts(dropna=False).items():
    lines.append(f"  {str(method):<40} {count:,}")
lines.append("")

lines.append("--- Skew Assessability ---")
lines.append(f"skew_assessable = True:         {merged['skew_assessable'].sum():,}")
lines.append(f"skew_assessable = False:        {(~merged['skew_assessable']).sum():,}")
lines.append("")

lines.append("--- Not Assessable Reasons ---")
for reason, count in merged["skew_not_assessable_reason"].value_counts(dropna=False).items():
    lines.append(f"  {str(reason):<55} {count:,}")
lines.append("")

lines.append("--- Skew Flags (assessable rows only) ---")
ass = merged[merged["skew_assessable"]]
lines.append(f"  3x flag:  {ass['skew_flag_3x'].sum():,}  ({100*ass['skew_flag_3x'].mean():.1f}% of assessable)")
lines.append(f"  5x flag:  {ass['skew_flag_5x'].sum():,}  ({100*ass['skew_flag_5x'].mean():.1f}% of assessable)")
lines.append(f"  10x flag: {ass['skew_flag_10x'].sum():,}  ({100*ass['skew_flag_10x'].mean():.1f}% of assessable)")
lines.append("")

lines.append("--- Zero min_AF flags (where max_AF >= 1%) ---")
lines.append(f"  zero_min_af_flag = True: {merged['zero_min_af_flag'].sum():,}")
lines.append("")

lines.append("--- AF Range Check (all values should be 0-1) ---")
for pop in PRIMARY_POPS:
    col = f"{pop}_af"
    bad = merged[(merged[col] < 0) | (merged[col] > 1)]
    lines.append(f"  {pop}: bad AF values = {len(bad)}")
lines.append("")

lines.append("--- Top max_pop (population with highest AF) ---")
lines.append(merged["max_pop"].value_counts(dropna=False).to_string())
lines.append("")
lines.append("--- Top min_pop (population with lowest AF) ---")
lines.append(merged["min_pop"].value_counts(dropna=False).to_string())
lines.append("")

lines.append("--- disparity_ratio distribution (assessable, non-zero-min) ---")
ratio_vals = merged[merged["skew_assessable"] & ~merged["zero_min_af_flag"]]["disparity_ratio"].dropna()
lines.append(f"  n values:   {len(ratio_vals):,}")
lines.append(f"  min ratio:  {ratio_vals.min():.3f}")
lines.append(f"  median:     {ratio_vals.median():.3f}")
lines.append(f"  95th pct:   {ratio_vals.quantile(0.95):.3f}")
lines.append(f"  max ratio:  {ratio_vals.max():.3f}")
lines.append("")

lines.append("--- Top 10 highest-ratio variants (spot-check these) ---")
top_ratio = merged[merged["skew_assessable"]].nlargest(10, "disparity_ratio")[
    ["rsID", "Gene", "Drug(s)", "Phenotype Category",
     "max_pop", "max_AF", "min_pop", "min_AF", "disparity_ratio"]
]
lines.append(top_ratio.to_string(index=False))
lines.append("")

lines.append("--- ClinVar star distribution among matched rows ---")
lines.append(merged[merged["in_clinvar"]]["clinvar_max_stars"].value_counts().sort_index().to_string())
lines.append("")

lines.append("--- pharmgkb_evidence_level distribution ---")
lines.append(merged["pharmgkb_evidence_level"].value_counts(dropna=False).to_string())

report = "\n".join(lines)
print("\n" + report)

# Save sanity check report
sanity_path = CLEAN / "sanity_checks.txt"
with open(sanity_path, "w") as f:
    f.write(report)
print(f"\nSanity report saved: {sanity_path}")

# ── Step 8: Save merged table ────────────────────────────────────────────────
out_path = CLEAN / "merged_table.csv"
merged.to_csv(out_path, index=False)
print(f"Saved merged table: {out_path}")
print(f"  Shape: {merged.shape}")
print("\nDone.")
