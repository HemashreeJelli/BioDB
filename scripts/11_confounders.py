"""
11_confounders.py -- Step 8: confounder controls for the ancestry analysis

Step 10 found that 45.4% of high-skew / low-evidence rsIDs have an evidence
base excluding the population where the allele peaks, vs 3.7% of high-skew /
high-evidence rsIDs. That contrast is confounded: high-evidence variants carry
a median of 32 published cohorts, low-evidence ones a median of 3, and a
variant with 3 cohorts is unlikely to cover any particular population by
chance. This step tests whether anything survives that confound.

Analyses:
  1. Mismatch rate as a function of cohort count (the raw confound).
  2. Permutation null: expected mismatch if each rsID's cohort ancestries were
     drawn at random from the global pool, holding its cohort count fixed.
  3. Logistic regression of mismatch on skew, evidence tier, cohort count,
     participant count, and curation recency.
  4. Participant-weighted definition of "includes the peak population".
  5. Missingness: is ancestry coverage itself related to skew and evidence?
  6. Drug-granularity overlap: rsIDs classified into both Group 1 and Group 2.
  7. Group X exclusion asymmetry in the baseline Fisher test.

Inputs:
  data_clean/final_table.csv
  results/robustness/evidence_ancestry.csv   (from step 10)
  data_clean/annotation_resolution_v3.csv
  variantAnnotations/study_parameters.tsv

Outputs (results/robustness/):
  confound_cohort_count.csv
  permutation_null.csv
  logistic_mismatch.csv
  weighted_ancestry_match.csv
  ancestry_missingness.csv
  group_overlap_rsids.csv
  groupx_exclusion.csv
  confounders_report.txt

No new dependencies: the logistic regression is fitted by Newton-Raphson and
Wald p-values come from the normal CDF via math.erfc.
"""

import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT   = Path(__file__).resolve().parent.parent
CLEAN  = ROOT / "data_clean"
VARANN = ROOT / "variantAnnotations"
OUT    = ROOT / "results" / "robustness"
OUT.mkdir(parents=True, exist_ok=True)

PRIMARY_POPS = ["AFR", "AMR", "EAS", "NFE", "SAS"]
N_PERM       = 1_000
SEED         = 20261025
WEIGHT_SHARE = 0.10   # peak population must hold >=10% of participants

ANCESTRY_MAP = {
    "European": "NFE",
    "East Asian": "EAS",
    "African American/Afro-Caribbean": "AFR",
    "Sub-Saharan African": "AFR",
    "Central/South Asian": "SAS",
    "Latino": "AMR",
    "American": "AMR",
}

lines = []


def say(text=""):
    print(text)
    lines.append(str(text))


# ── Logistic regression (Newton-Raphson, no statsmodels) ────────────────────
def norm_sf_two_sided(z):
    return math.erfc(abs(z) / math.sqrt(2.0))


def logistic_fit(X, y, names, ridge=1e-6, max_iter=100, tol=1e-10):
    """Fit logistic regression. Returns a DataFrame of coefficients."""
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    n, k = X.shape
    beta = np.zeros(k)
    for _ in range(max_iter):
        eta = np.clip(X @ beta, -35, 35)
        p = 1.0 / (1.0 + np.exp(-eta))
        W = np.clip(p * (1 - p), 1e-10, None)
        grad = X.T @ (y - p) - ridge * beta
        H = (X * W[:, None]).T @ X + ridge * np.eye(k)
        try:
            step = np.linalg.solve(H, grad)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(H, grad, rcond=None)[0]
        beta_new = beta + step
        if np.max(np.abs(beta_new - beta)) < tol:
            beta = beta_new
            break
        beta = beta_new
    eta = np.clip(X @ beta, -35, 35)
    p = 1.0 / (1.0 + np.exp(-eta))
    W = np.clip(p * (1 - p), 1e-10, None)
    H = (X * W[:, None]).T @ X + ridge * np.eye(k)
    cov = np.linalg.pinv(H)
    se = np.sqrt(np.maximum(np.diag(cov), 0))
    z = np.divide(beta, se, out=np.zeros_like(beta), where=se > 0)
    return pd.DataFrame({
        "term": names,
        "coef": beta.round(4),
        "std_err": se.round(4),
        "odds_ratio": np.exp(beta).round(3),
        "or_ci_low": np.exp(beta - 1.96 * se).round(3),
        "or_ci_high": np.exp(beta + 1.96 * se).round(3),
        "z": z.round(3),
        "p_value": [norm_sf_two_sided(v) for v in z],
        "n_obs": n,
    })


# ── Load ────────────────────────────────────────────────────────────────────
say("=== CONFOUNDER CONTROLS FOR THE ANCESTRY ANALYSIS ===")
say()

df = pd.read_csv(CLEAN / "final_table.csv", low_memory=False)
ea = pd.read_csv(OUT / "evidence_ancestry.csv")
sp = pd.read_csv(VARANN / "study_parameters.tsv", sep="\t", low_memory=False)
res = pd.read_csv(CLEAN / "annotation_resolution_v3.csv", low_memory=False)

sp["aid"] = pd.to_numeric(sp["Variant Annotation ID"], errors="coerce")
res["aid"] = pd.to_numeric(res["annotation_id"], errors="coerce")
sp["group_label"] = (sp["Biogeographical Groups"].fillna("")
                     .astype(str).str.split(",").str[0].str.strip())
sp["mapped_pop"] = sp["group_label"].map(ANCESTRY_MAP)
sp["participants"] = (pd.to_numeric(sp["Study Cases"], errors="coerce").fillna(0)
                      + pd.to_numeric(sp["Study Controls"], errors="coerce").fillna(0))

studies = sp.merge(res[["aid", "rsID"]].drop_duplicates(), on="aid", how="inner")
say(f"studies joined: {len(studies):,}   rsIDs: {studies['rsID'].nunique():,}")

ea["excludes_peak"] = (ea["ancestry_match"] == "evidence_excludes_peak_population").astype(int)
usable = ea[ea["has_ancestry_data"]].copy()
usable["stratum"] = (np.where(usable["high_skew"], "high skew", "low skew")
                     + " / " + usable["tier"] + " evidence")
say(f"rsIDs with usable ancestry data: {len(usable):,}")
say()


# ── 1. The raw confound ─────────────────────────────────────────────────────
say("--- 1. Mismatch rate by number of published cohorts ---")
bins = [0, 1, 2, 5, 10, 30, 10**6]
labels = ["1", "2", "3-5", "6-10", "11-30", ">30"]
usable["cohort_bin"] = pd.cut(usable["n_studies"], bins=bins, labels=labels)
conf = (usable.groupby("cohort_bin", observed=True)
        .agg(rsids=("excludes_peak", "size"),
             pct_excluding_peak=("excludes_peak", lambda s: round(100 * s.mean(), 1)))
        .reset_index())
say(conf.to_string(index=False))
conf.to_csv(OUT / "confound_cohort_count.csv", index=False)
say()
say("  Median cohorts per stratum:")
med = usable.groupby("stratum")["n_studies"].median()
for k, v in med.items():
    say(f"    {k:<34} {v:>6.0f}")
say()
say("  Mismatch falls monotonically as cohort count rises, and the two strata")
say("  differ ~10-fold in cohort count, so the headline contrast is confounded.")
say()


# ── 2. Permutation null ─────────────────────────────────────────────────────
# For each rsID, hold the number of distinct study populations fixed and draw
# that many populations from the global marginal distribution of study
# populations. Mismatch under this null is what cohort scarcity alone predicts.
say(f"--- 2. Permutation null ({N_PERM:,} draws) ---")

pop_sets = (studies.dropna(subset=["mapped_pop"])
            .groupby("rsID")["mapped_pop"].apply(lambda s: sorted(set(s))))
global_counts = (studies.dropna(subset=["mapped_pop"])
                 .groupby("rsID")["mapped_pop"].apply(lambda s: sorted(set(s)))
                 .explode().value_counts())
pops = list(global_counts.index)
weights = (global_counts / global_counts.sum()).reindex(pops).to_numpy()
say("  Global cohort-population marginals (share of rsIDs covering each):")
for p, w in zip(pops, weights):
    say(f"    {p}: {100*w:.1f}%")

rng = np.random.default_rng(SEED)
sub = usable[usable["rsID"].isin(pop_sets.index)].copy()
sub["k_pops"] = sub["rsID"].map(pop_sets.apply(len))

obs_rate = {}
exp_rate = {}
for stratum, grp in sub.groupby("stratum"):
    peaks = grp["max_pop"].to_numpy()
    ks = grp["k_pops"].to_numpy()
    hits = 0
    for _ in range(N_PERM):
        miss = 0
        for peak, k in zip(peaks, ks):
            k_eff = min(int(k), len(pops))
            drawn = rng.choice(pops, size=k_eff, replace=False,
                               p=weights / weights.sum())
            if peak not in drawn:
                miss += 1
        hits += miss / len(grp)
    exp_rate[stratum] = 100 * hits / N_PERM
    obs_rate[stratum] = 100 * grp["excludes_peak"].mean()

perm = pd.DataFrame({
    "stratum": list(obs_rate.keys()),
    "n_rsids": [int((sub["stratum"] == s).sum()) for s in obs_rate],
    "observed_pct_excluding": [round(obs_rate[s], 1) for s in obs_rate],
    "null_expected_pct": [round(exp_rate[s], 1) for s in obs_rate],
})
perm["observed_minus_null"] = (perm["observed_pct_excluding"]
                               - perm["null_expected_pct"]).round(1)
say()
say(perm.to_string(index=False))
perm.to_csv(OUT / "permutation_null.csv", index=False)
say()
say("  A difference near zero means cohort scarcity and the overall skew of")
say("  the literature already explain the mismatch, with no ancestry-specific")
say("  effect left over.")
say()


# ── 3. Logistic regression ──────────────────────────────────────────────────
say("--- 3. Logistic regression of mismatch ---")

part = studies.groupby("rsID")["participants"].sum().rename("total_participants")
recency = (df.groupby("rsID")["Latest History Date (YYYY-MM-DD)"]
           .max().rename("latest_date"))

reg = usable.merge(part, on="rsID", how="left").merge(recency, on="rsID", how="left")
reg["total_participants"] = reg["total_participants"].fillna(0)
reg["curation_year"] = pd.to_datetime(reg["latest_date"], errors="coerce").dt.year
reg["curation_year"] = reg["curation_year"].fillna(reg["curation_year"].median())

reg = reg.dropna(subset=["n_studies", "max_pop"]).copy()
reg["log_cohorts"] = np.log1p(reg["n_studies"])
reg["log_participants"] = np.log1p(reg["total_participants"])
reg["high_skew_i"] = reg["high_skew"].astype(int)
reg["high_ev_i"] = (reg["tier"] == "high").astype(int)
reg["year_c"] = reg["curation_year"] - reg["curation_year"].mean()

specs = [
    ("skew + evidence only",
     ["high_skew_i", "high_ev_i"]),
    ("+ cohort count",
     ["high_skew_i", "high_ev_i", "log_cohorts"]),
    ("full model",
     ["high_skew_i", "high_ev_i", "log_cohorts", "log_participants", "year_c"]),
]

all_fits = []
for name, cols in specs:
    X = np.column_stack([np.ones(len(reg))] + [reg[c].to_numpy(float) for c in cols])
    fit = logistic_fit(X, reg["excludes_peak"].to_numpy(float), ["intercept"] + cols)
    fit.insert(0, "model", name)
    all_fits.append(fit)
    say(f"  Model: {name}  (n = {len(reg):,})")
    for _, r in fit.iterrows():
        star = "  *" if r["p_value"] < 0.05 else ""
        say(f"    {r['term']:<20} OR = {r['odds_ratio']:>8.3f}  "
            f"95% CI [{r['or_ci_low']:.3f}, {r['or_ci_high']:.3f}]   "
            f"p = {r['p_value']:.3g}{star}")
    say()

say(f"  High-evidence rsIDs in the model: "
    f"{int((reg['tier'] == 'high').sum())} of {len(reg):,}. The evidence-tier")
say("  coefficient is therefore poorly powered: a wide CI spanning 1 means")
say("  'no detectable effect', not 'no effect'.")
say()

pd.concat(all_fits, ignore_index=True).to_csv(OUT / "logistic_mismatch.csv", index=False)
say("  Read the coefficient on high_skew_i in the full model: if it loses")
say("  significance once cohort count is included, the ancestry contrast is")
say("  explained by how much literature exists, not by skew.")
say()


# ── 4. Participant-weighted "includes peak" ─────────────────────────────────
say(f"--- 4. Participant-weighted match (peak needs >={100*WEIGHT_SHARE:.0f}% of participants) ---")
pw = (studies.dropna(subset=["mapped_pop"])
      .groupby(["rsID", "mapped_pop"])["participants"].sum().reset_index())
tot = pw.groupby("rsID")["participants"].sum().rename("tot")
pw = pw.merge(tot, on="rsID")
pw["share"] = np.where(pw["tot"] > 0, pw["participants"] / pw["tot"], np.nan)

peak_share = (usable[["rsID", "max_pop", "stratum", "excludes_peak"]]
              .merge(pw, left_on=["rsID", "max_pop"],
                     right_on=["rsID", "mapped_pop"], how="left"))
peak_share["share"] = peak_share["share"].fillna(0.0)
peak_share["weighted_excludes"] = (peak_share["share"] < WEIGHT_SHARE).astype(int)

wsum = (peak_share.groupby("stratum")
        .agg(rsids=("weighted_excludes", "size"),
             pct_excl_binary=("excludes_peak", lambda s: round(100 * s.mean(), 1)),
             pct_excl_weighted=("weighted_excludes", lambda s: round(100 * s.mean(), 1)))
        .reset_index())
say(wsum.to_string(index=False))
wsum.to_csv(OUT / "weighted_ancestry_match.csv", index=False)
say()
say("  NOTE: participant counts are missing or zero for many annotations, so")
say("  the weighted figure is an upper bound on exclusion.")
say()


# ── 5. Missingness of ancestry data ─────────────────────────────────────────
say("--- 5. Is ancestry coverage itself related to skew and evidence? ---")
peak_all = (df[df["skew_assessable"]].groupby("rsID")
            .agg(high_skew=("skew_flag_5x", "any"),
                 tier=("pgkb_evidence_tier",
                       lambda s: "high" if (s == "high").any() else "low"))
            .reset_index())
miss = peak_all.merge(ea[["rsID", "has_ancestry_data"]], on="rsID", how="left")
miss["has_ancestry_data"] = miss["has_ancestry_data"].fillna(False)
miss["stratum"] = (np.where(miss["high_skew"], "high skew", "low skew")
                   + " / " + miss["tier"] + " evidence")
msum = (miss.groupby("stratum")
        .agg(rsids=("has_ancestry_data", "size"),
             with_ancestry=("has_ancestry_data", "sum"),
             pct_covered=("has_ancestry_data", lambda s: round(100 * s.mean(), 1)))
        .reset_index())
say(msum.to_string(index=False))
msum.to_csv(OUT / "ancestry_missingness.csv", index=False)
say()
say("  Coverage is NOT missing at random: the stratum carrying the headline")
say("  claim is the least well covered, so its mismatch rate rests on a")
say("  non-random subset.")
say()


# ── 6. Drug-granularity overlap ─────────────────────────────────────────────
say("--- 6. rsIDs classified into BOTH Group 1 and Group 2 ---")
assess = df[df["skew_assessable"]]
by_rsid = assess.groupby("rsID")["group_pgkb"].apply(set)
both = [r for r, s in by_rsid.items()
        if any("Group 1" in x for x in s) and any("Group 2" in x for x in s)]
ov = assess[assess["rsID"].isin(both)][
    ["rsID", "Gene", "Drug(s)", "pharmgkb_evidence_level", "group_pgkb",
     "max_AF", "min_AF", "absolute_AF_diff", "disparity_ratio"]
].sort_values(["rsID", "pharmgkb_evidence_level"])
ov.to_csv(OUT / "group_overlap_rsids.csv", index=False)

g2_rows = assess[assess["group_pgkb"].str.startswith("Group 2")]
g2_in_both = g2_rows[g2_rows["rsID"].isin(both)]
say(f"  rsIDs in both Group 1 and Group 2: {len(both)}")
say(f"    {', '.join(sorted(both))}")
say(f"  Group 2 rows belonging to those rsIDs: {len(g2_in_both):,} "
    f"({100*len(g2_in_both)/len(g2_rows):.1f}% of Group 2)")
say(f"  Group 2 rsIDs that are high-evidence elsewhere: {len(both)} of "
    f"{g2_rows['rsID'].nunique()} ({100*len(both)/g2_rows['rsID'].nunique():.1f}%)")
say()
say("  For these variants the frequency data are identical across groups and")
say("  only the drug differs, so the gap is drug-phenotype curation rather")
say("  than a variant nobody has studied.")
say()


# ── 7. Group X exclusion asymmetry ──────────────────────────────────────────
say("--- 7. Group X exclusion asymmetry in the baseline test ---")
gx = []
for tier in ["high", "low"]:
    t = df[df["pgkb_evidence_tier"] == tier]
    gx.append({
        "evidence_tier": tier,
        "rows_total": len(t),
        "rows_assessable": int(t["skew_assessable"].sum()),
        "rows_group_x": int((~t["skew_assessable"]).sum()),
        "pct_excluded": round(100 * (~t["skew_assessable"]).mean(), 1),
    })
gxdf = pd.DataFrame(gx)
say(gxdf.to_string(index=False))
gxdf.to_csv(OUT / "groupx_exclusion.csv", index=False)
say()
say("  The baseline Fisher test therefore compares the common, assessable")
say("  minority of high-evidence rows against nearly all low-evidence rows.")
say()


report = OUT / "confounders_report.txt"
report.write_text("\n".join(lines), encoding="utf-8")
print(f"\nSaved report: {report}")
print(f"Saved tables to: {OUT}")
print("\nDone.")
