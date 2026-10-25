"""
10_robustness.py -- Step 7: robustness, baseline, and evidence-ancestry analyses

Addresses reviewer criticisms of the four-group audit. Nothing here changes the
pipeline's outputs; it recomputes the headline numbers under alternative
assumptions and adds one new analysis (evidence ancestry).

Inputs:
  data_clean/final_table.csv               -- classified table from step 09
  data_clean/annotation_resolution_v3.csv  -- annotation_id -> rsID link
  variantAnnotations/study_parameters.tsv  -- per-study Biogeographical Groups

Outputs (results/robustness/):
  baseline_association.csv     -- skew rate WITHIN each evidence tier
  delta_af_by_group.csv        -- absolute AF difference per group
  ac_fragility.csv             -- allele counts behind high-skew calls
  detection_floor.csv          -- rows whose min AF sits near 1/AN
  sensitivity_arms.csv         -- Group 2 share under each robustness arm
  downsample_arm.csv           -- skew recomputed at a common allele number
  rsid_level_groups.csv        -- groups counted once per rsID
  skew_direction.csv           -- which population carries the HIGH frequency
  evidence_ancestry.csv        -- per-rsID study ancestry vs. frequency peak
  ancestry_mismatch_summary.csv
  robustness_report.txt        -- everything above, printed

Design notes:
  - Fisher's exact test is implemented locally (math.lgamma) so this script
    adds no dependency on scipy.
  - Wilson score intervals are used for allele frequency CIs; the
    "conservative ratio" divides the lower bound of max_AF by the upper
    bound of min_AF, so it cannot be inflated by sampling noise alone.
  - Downsampling draws binomial counts at a common AN with a fixed seed,
    which tests whether extreme ratios survive equal sample sizes.
"""

import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT    = Path(__file__).resolve().parent.parent
CLEAN   = ROOT / "data_clean"
VARANN  = ROOT / "variantAnnotations"
OUT     = ROOT / "results" / "robustness"
OUT.mkdir(parents=True, exist_ok=True)

PRIMARY_POPS = ["AFR", "AMR", "EAS", "NFE", "SAS"]
MIN_AN       = 1_000
MIN_MAX_AF   = 0.01
SKEW_THRESH  = 5.0
MIN_AC_FILTER = 10      # sensitivity arm: require >10 copies in the low population
DOWNSAMPLE_AN = 30_000  # common allele number for the equal-sample-size arm
N_REPLICATES  = 200
SEED          = 20261025

# PharmGKB biogeographical group -> gnomAD primary population.
# Near Eastern maps to MID and Oceanian has no gnomAD primary equivalent, so
# both are left unmapped rather than forced into one of the five groups.
ANCESTRY_MAP = {
    "European": "NFE",
    "East Asian": "EAS",
    "African American/Afro-Caribbean": "AFR",
    "Sub-Saharan African": "AFR",
    "Central/South Asian": "SAS",
    "Latino": "AMR",
    "American": "AMR",
}
UNMAPPED = {"Near Eastern", "Oceanian", "Custom", "Unknown", ""}

lines = []


def say(text=""):
    print(text)
    lines.append(str(text))


# ── Statistics helpers ──────────────────────────────────────────────────────
def _log_hypergeom(a, b, c, d):
    """log P of one 2x2 table under fixed margins."""
    n = a + b + c + d
    lg = math.lgamma
    return (lg(a + b + 1) + lg(c + d + 1) + lg(a + c + 1) + lg(b + d + 1)
            - lg(a + 1) - lg(b + 1) - lg(c + 1) - lg(d + 1) - lg(n + 1))


def fisher_exact_2x2(a, b, c, d):
    """Two-sided Fisher's exact test. Returns (odds_ratio, p_value)."""
    a, b, c, d = int(a), int(b), int(c), int(d)
    observed = math.exp(_log_hypergeom(a, b, c, d))
    row1, row2 = a + b, c + d
    col1 = a + c
    n = row1 + row2
    lo = max(0, col1 - row2)
    hi = min(row1, col1)
    p = 0.0
    for x in range(lo, hi + 1):
        y = row1 - x
        z = col1 - x
        w = row2 - z
        prob = math.exp(_log_hypergeom(x, y, z, w))
        if prob <= observed * (1 + 1e-9):
            p += prob
    odds = float("inf") if (b == 0 or c == 0) else (a * d) / (b * c)
    return odds, min(p, 1.0)


def wilson(count, total, z=1.96):
    """Wilson score interval for a proportion. Robust at count = 0."""
    if not total or total <= 0 or pd.isna(count) or pd.isna(total):
        return (np.nan, np.nan)
    p = count / total
    denom = 1 + z**2 / total
    centre = (p + z**2 / (2 * total)) / denom
    half = (z / denom) * math.sqrt(p * (1 - p) / total + z**2 / (4 * total**2))
    return (max(0.0, centre - half), min(1.0, centre + half))


# ── Load ────────────────────────────────────────────────────────────────────
say("=== ROBUSTNESS & BASELINE ANALYSES ===")
say()

df = pd.read_csv(CLEAN / "final_table.csv", low_memory=False)
say(f"Loaded final_table: {len(df):,} rows, {df['rsID'].nunique():,} unique rsIDs")

assessable = df[df["skew_assessable"]].copy()
say(f"Assessable rows: {len(assessable):,}")
say()


# ── 1. Baseline: skew rate WITHIN each evidence tier ────────────────────────
# The headline "30.7% of rows are high-skew / low-evidence" is driven largely
# by the fact that 95.2% of all rows are Level 3 or 4. The informative
# quantity is the skew RATE inside each tier.
say("--- 1. Baseline association: skew rate within evidence tier ---")

rows = []
for level, label in [("row", "rsID x gene x drug rows"), ("rsid", "unique rsIDs")]:
    if level == "row":
        tab = pd.crosstab(assessable["pgkb_evidence_tier"], assessable["skew_flag_5x"])
    else:
        per_rsid = assessable.groupby("rsID").agg(
            tier=("pgkb_evidence_tier", lambda s: "high" if (s == "high").any() else "low"),
            high_skew=("skew_flag_5x", "any"),
        )
        tab = pd.crosstab(per_rsid["tier"], per_rsid["high_skew"])

    tab = tab.reindex(index=["high", "low"], columns=[False, True], fill_value=0)
    a, b = tab.loc["high", True], tab.loc["high", False]
    c, d = tab.loc["low", True], tab.loc["low", False]
    odds, p = fisher_exact_2x2(a, b, c, d)

    hi_rate = 100 * a / (a + b) if (a + b) else np.nan
    lo_rate = 100 * c / (c + d) if (c + d) else np.nan
    say(f"  Unit: {label}")
    say(f"    High evidence: {a:>5} / {a+b:>5} high-skew  ({hi_rate:.1f}%)")
    say(f"    Low evidence:  {c:>5} / {c+d:>5} high-skew  ({lo_rate:.1f}%)")
    say(f"    Fisher's exact: OR = {odds:.2f}, p = {p:.3g}")
    rows.append({
        "unit": label,
        "high_ev_highskew": a, "high_ev_total": a + b, "high_ev_pct": round(hi_rate, 1),
        "low_ev_highskew": c, "low_ev_total": c + d, "low_ev_pct": round(lo_rate, 1),
        "odds_ratio": round(odds, 3), "p_value": p,
    })

pd.DataFrame(rows).to_csv(OUT / "baseline_association.csv", index=False)
say("  NOTE: high-evidence variants are MORE often skewed, so a raw count of")
say("        high-skew/low-evidence rows cannot by itself show neglect.")
say()


# ── 2. Absolute AF difference by group ──────────────────────────────────────
# Guards against the criticism that large ratios are artifacts of tiny
# denominators: absolute differences do not have that failure mode.
say("--- 2. Absolute AF difference by group ---")
dstats = (assessable.groupby("group_pgkb")["absolute_AF_diff"]
          .describe()[["count", "25%", "50%", "75%", "max"]].round(4))
say(dstats.to_string())
dstats.to_csv(OUT / "delta_af_by_group.csv")
say()


# ── 3. Fragility of the ratio: allele counts behind high-skew calls ─────────
say("--- 3. Allele counts in the low-frequency population (high-skew rows) ---")
high_skew = assessable[assessable["skew_flag_5x"]].copy()


def low_pop_field(row, suffix):
    pop = row["min_pop"]
    return row.get(f"{pop}_{suffix}", np.nan) if isinstance(pop, str) else np.nan


high_skew["min_pop_ac"] = high_skew.apply(low_pop_field, suffix="ac", axis=1)
high_skew["min_pop_an"] = high_skew.apply(low_pop_field, suffix="an", axis=1)

frag = []
for thr in [0, 1, 5, 10, 20, 50, 100]:
    n = int((high_skew["min_pop_ac"] > thr).sum())
    frag.append({"min_ac_threshold": thr, "high_skew_rows_retained": n,
                 "pct_retained": round(100 * n / len(high_skew), 1)})
    say(f"  AC > {thr:>3}: {n:>5} of {len(high_skew):,} retained "
        f"({100*n/len(high_skew):.1f}%)")
say(f"  median AC in low population: {high_skew['min_pop_ac'].median():,.0f}")
say(f"  ratio > 1000 AND AC <= 10:   "
    f"{int(((high_skew['disparity_ratio'] > 1000) & (high_skew['min_pop_ac'] <= 10)).sum())}")
pd.DataFrame(frag).to_csv(OUT / "ac_fragility.csv", index=False)
say()


# ── 4. Detection floor and confidence-interval-aware ratios ─────────────────
# A population with a small AN cannot resolve very low frequencies. Rows whose
# min_AF sits within a few multiples of 1/AN are measuring the floor, not biology.
say("--- 4. Detection floor and conservative (CI-based) ratios ---")
high_skew["min_pop_detect_floor"] = 1.0 / high_skew["min_pop_an"]
high_skew["floor_ratio"] = high_skew["min_AF"] / high_skew["min_pop_detect_floor"]
near_floor = high_skew["floor_ratio"] < 5
say(f"  high-skew rows whose min AF is < 5x its detection floor: "
    f"{int(near_floor.sum())} ({100*near_floor.mean():.1f}%)")


def conservative_ratio(row):
    """Lower bound of max_AF divided by upper bound of min_AF."""
    hi_pop, lo_pop = row["max_pop"], row["min_pop"]
    if not isinstance(hi_pop, str) or not isinstance(lo_pop, str):
        return np.nan
    hi_lo, _ = wilson(row.get(f"{hi_pop}_ac"), row.get(f"{hi_pop}_an"))
    _, lo_hi = wilson(row.get(f"{lo_pop}_ac"), row.get(f"{lo_pop}_an"))
    if pd.isna(hi_lo) or pd.isna(lo_hi) or lo_hi <= 0:
        return np.nan
    return hi_lo / lo_hi


high_skew["conservative_ratio"] = high_skew.apply(conservative_ratio, axis=1)
survives = (high_skew["conservative_ratio"] >= SKEW_THRESH).sum()
say(f"  high-skew rows still >= {SKEW_THRESH:g}x using CI bounds: "
    f"{int(survives)} of {len(high_skew):,} ({100*survives/len(high_skew):.1f}%)")
say(f"  median conservative ratio: {high_skew['conservative_ratio'].median():.1f} "
    f"(median point estimate: {high_skew['disparity_ratio'].median():.1f})")

high_skew[["rsID", "Gene", "Drug(s)", "max_pop", "max_AF", "min_pop", "min_AF",
           "min_pop_ac", "min_pop_an", "disparity_ratio", "conservative_ratio",
           "floor_ratio"]].to_csv(OUT / "detection_floor.csv", index=False)
say()


# ── 5. Sensitivity arms: does the headline share survive? ───────────────────
say("--- 5. Group 2 share under each robustness arm ---")
arms = []


def group2_share(frame, name, note=""):
    n = len(frame)
    g2 = int(frame["group_pgkb"].astype(str).str.startswith("Group 2").sum())
    pct = round(100 * g2 / n, 1) if n else np.nan
    arms.append({"arm": name, "rows": n, "group2_rows": g2,
                 "group2_pct": pct, "note": note})
    say(f"  {name:<42} {g2:>5} / {n:<6} ({pct:.1f}%)  {note}")


group2_share(df, "primary (as published)")
group2_share(df[df["allele_selection_method"] != "highest_ac_fallback"],
             "excluding highest-AC fallback rows", "allele assignment check")

keep_ac = set(high_skew.loc[high_skew["min_pop_ac"] > MIN_AC_FILTER].index)
ac_arm = df.loc[[i for i in df.index
                 if i in keep_ac or i not in high_skew.index]]
group2_share(ac_arm, f"requiring AC > {MIN_AC_FILTER} in low population",
             "ratio fragility check")

ci_keep = set(high_skew.loc[high_skew["conservative_ratio"] >= SKEW_THRESH].index)
ci_arm = df.loc[[i for i in df.index
                 if i in ci_keep or i not in high_skew.index]]
group2_share(ci_arm, "requiring CI-based ratio >= 5x", "sampling-noise check")

pd.DataFrame(arms).to_csv(OUT / "sensitivity_arms.csv", index=False)
say()


# ── 6. Equal-sample-size (downsampling) arm ─────────────────────────────────
# NFE allele numbers are an order of magnitude larger than EAS, so the most
# extreme ratios are easier to produce when NFE is the low population. This
# arm redraws every population at a common AN and re-tests the 5x threshold.
say(f"--- 6. Downsampling to a common AN of {DOWNSAMPLE_AN:,} "
    f"({N_REPLICATES} replicates) ---")

rng = np.random.default_rng(SEED)
sub = assessable.reset_index(drop=True)
af_matrix = {p: sub[f"{p}_af"].to_numpy(dtype=float) for p in PRIMARY_POPS}
an_matrix = {p: sub[f"{p}_an"].to_numpy(dtype=float) for p in PRIMARY_POPS}
eligible = {p: (~np.isnan(af_matrix[p])) & (an_matrix[p] >= MIN_AN) for p in PRIMARY_POPS}

flag_counts = np.zeros(len(sub))
for _ in range(N_REPLICATES):
    drawn = {}
    for p in PRIMARY_POPS:
        af = np.nan_to_num(af_matrix[p], nan=0.0)
        counts = rng.binomial(DOWNSAMPLE_AN, np.clip(af, 0, 1))
        vals = counts / DOWNSAMPLE_AN
        vals[~eligible[p]] = np.nan
        drawn[p] = vals
    stacked = np.vstack([drawn[p] for p in PRIMARY_POPS])
    with np.errstate(invalid="ignore"):
        mx = np.nanmax(stacked, axis=0)
        mn = np.nanmin(stacked, axis=0)
    # A zero that only appears after downsampling is a censored measurement,
    # not evidence of absence, so it does NOT count as high skew here. The
    # full-data pipeline treats a true zero as high skew, but at AN = 30,000 a
    # frequency of 2e-5 rounds to zero copies most of the time, which would
    # manufacture high-skew calls out of reduced sampling depth.
    flagged = np.where(
        (mx < MIN_MAX_AF) | (mn == 0), False,
        mx / np.where(mn == 0, np.nan, mn) >= SKEW_THRESH,
    )
    flag_counts += np.nan_to_num(flagged.astype(float))

sub["downsample_highskew_frac"] = flag_counts / N_REPLICATES
orig_high = sub["skew_flag_5x"].astype(bool)
stable = sub["downsample_highskew_frac"] >= 0.95
say(f"  originally high-skew rows:                 {int(orig_high.sum()):,}")
say(f"  still high-skew in >=95% of replicates:    "
    f"{int((orig_high & stable).sum()):,}")
say(f"  lost under equal sample sizes:             "
    f"{int((orig_high & ~stable).sum()):,}")
say(f"  newly high-skew (was low):                 "
    f"{int((~orig_high & stable).sum()):,}")
say(f"  preservation rate:                         "
    f"{100*(orig_high & stable).sum()/max(int(orig_high.sum()),1):.1f}%")
say("  Downsample-induced zeros are treated as censored, not as absence; the")
say("  arm therefore measures only skew that survives equal sampling depth.")
sub[["rsID", "Gene", "Drug(s)", "group_pgkb", "max_pop", "min_pop",
     "disparity_ratio", "skew_flag_5x", "downsample_highskew_frac"]].to_csv(
    OUT / "downsample_arm.csv", index=False)
say()


# ── 7. rsID-level group counts ──────────────────────────────────────────────
# Rows are not independent: 3,834 rows come from 2,506 rsIDs.
say("--- 7. Groups counted once per rsID ---")
per_rsid = df.groupby("rsID").agg(
    tier=("pgkb_evidence_tier", lambda s: "high" if (s == "high").any() else "low"),
    skew=("skew_tier", lambda s: "high" if (s == "high").any()
          else ("low" if (s == "low").any() else "not_assessable")),
)
rl = (per_rsid.groupby(["skew", "tier"]).size()
      .rename("unique_rsids").reset_index())
rl["pct_of_rsids"] = (100 * rl["unique_rsids"] / len(per_rsid)).round(1)
say(rl.to_string(index=False))
rl.to_csv(OUT / "rsid_level_groups.csv", index=False)
say(f"  total rsIDs: {len(per_rsid):,}")
say()


# ── 8. Direction of skew ────────────────────────────────────────────────────
# Which population carries the HIGH frequency matters for interpretation: the
# minimum-AF population is where the allele is RARE, so it is not the group
# most affected by that variant.
say("--- 8. Direction of skew: population carrying the HIGH frequency ---")
direction = (pd.crosstab(df["group_pgkb"], df["max_pop"])
             .reindex(columns=PRIMARY_POPS, fill_value=0))
say(direction.to_string())
direction.to_csv(OUT / "skew_direction.csv")
say()
say("  Median allele number per population (sampling asymmetry):")
for p in PRIMARY_POPS:
    say(f"    {p}: {df[f'{p}_an'].median():>10,.0f}")
say()


# ── 9. Evidence ancestry vs. frequency peak ─────────────────────────────────
# New analysis. Tests the premise motivating the whole audit: is the evidence
# behind a variant drawn from the population in which the allele is common?
say("--- 9. Evidence ancestry vs. population carrying the high frequency ---")

sp = pd.read_csv(VARANN / "study_parameters.tsv", sep="\t", low_memory=False)
res = pd.read_csv(CLEAN / "annotation_resolution_v3.csv", low_memory=False)

sp["aid"] = pd.to_numeric(sp["Variant Annotation ID"], errors="coerce")
res["aid"] = pd.to_numeric(res["annotation_id"], errors="coerce")

# "Multiple groups, <free text>" and similar carry a description after the
# first comma; keep only the leading label.
sp["group_label"] = (sp["Biogeographical Groups"].fillna("")
                     .astype(str).str.split(",").str[0].str.strip())
sp["mapped_pop"] = sp["group_label"].map(ANCESTRY_MAP)

studies = sp.merge(res[["aid", "rsID"]].drop_duplicates(), on="aid", how="inner")
say(f"  studies joined to an rsID: {len(studies):,} "
    f"({studies['rsID'].nunique():,} rsIDs)")


def summarise_ancestry(group):
    mapped = group["mapped_pop"].dropna()
    multi = (group["group_label"] == "Multiple groups").sum()
    unknown = group["group_label"].isin(UNMAPPED).sum()
    return pd.Series({
        "n_studies": len(group),
        "n_mapped_studies": len(mapped),
        "n_multiple_groups": int(multi),
        "n_unknown_or_unmapped": int(unknown),
        "study_pops": ";".join(sorted(mapped.unique())) if len(mapped) else "",
        "dominant_study_pop": mapped.value_counts().idxmax() if len(mapped) else np.nan,
    })


anc = studies.groupby("rsID").apply(summarise_ancestry, include_groups=False).reset_index()

peak = (df[df["skew_assessable"]]
        .groupby("rsID")
        .agg(max_pop=("max_pop", "first"),
             max_AF=("max_AF", "first"),
             high_skew=("skew_flag_5x", "any"),
             tier=("pgkb_evidence_tier",
                   lambda s: "high" if (s == "high").any() else "low"))
        .reset_index())

ea = peak.merge(anc, on="rsID", how="left")
ea["has_ancestry_data"] = ea["n_mapped_studies"].fillna(0) > 0


def mismatch(row):
    if not row["has_ancestry_data"]:
        return "no_ancestry_data"
    pops = set(str(row["study_pops"]).split(";")) - {""}
    if row["max_pop"] in pops:
        return "evidence_includes_peak_population"
    return "evidence_excludes_peak_population"


ea["ancestry_match"] = ea.apply(mismatch, axis=1)
ea.to_csv(OUT / "evidence_ancestry.csv", index=False)

say()
say("  Does the supporting evidence include the population where the allele peaks?")
overall = ea["ancestry_match"].value_counts()
for k, v in overall.items():
    say(f"    {k:<40} {v:>5} ({100*v/len(ea):.1f}%)")

say()
say("  Broken down by skew and evidence tier (rsIDs with ancestry data):")
sub_ea = ea[ea["has_ancestry_data"]].copy()
sub_ea["stratum"] = np.where(sub_ea["high_skew"], "high skew", "low skew") + \
                    " / " + sub_ea["tier"] + " evidence"
summary = (pd.crosstab(sub_ea["stratum"], sub_ea["ancestry_match"])
           .reindex(columns=["evidence_includes_peak_population",
                             "evidence_excludes_peak_population"],
                    fill_value=0))
summary["pct_excluding_peak"] = (
    100 * summary["evidence_excludes_peak_population"] / summary.sum(axis=1)
).round(1)
say(summary.to_string())
summary.to_csv(OUT / "ancestry_mismatch_summary.csv")

say()
say("  Peak population of rsIDs whose evidence EXCLUDES that population:")
excl = sub_ea[sub_ea["ancestry_match"] == "evidence_excludes_peak_population"]
say(excl["max_pop"].value_counts().to_string())

say()
say("  Dominant study population, high-skew rsIDs:")
say(sub_ea[sub_ea["high_skew"]]["dominant_study_pop"].value_counts().to_string())
say()

say("  CAVEAT: 'Multiple groups' and 'Unknown' labels carry no usable single")
say("  population, and PharmGKB's biogeographical labels are coarser than")
say("  gnomAD's genetic ancestry groups. Near Eastern and Oceanian studies are")
say("  left unmapped because neither has a gnomAD primary-population analogue.")
say()


# ── Save report ─────────────────────────────────────────────────────────────
report_path = OUT / "robustness_report.txt"
report_path.write_text("\n".join(lines), encoding="utf-8")
print(f"\nSaved report: {report_path}")
print(f"Saved tables to: {OUT}")
print("\nDone.")
