"""
12_coverage_and_clustering.py -- Step 9: continuous coverage metric and
                                 cluster-robust inference

Two weaknesses in the step 10/11 ancestry analysis are addressed here.

1. The binary "peak population excluded" metric is arbitrary when population
   frequencies are close. A variant at 28% in EAS and 27% in NFE is not
   meaningfully unstudied if only NFE cohorts exist. This step replaces the
   binary with a continuous measure:

       coverage_ratio = max(AF among studied populations) / peak AF

   and a threshold definition: a variant is covered if any studied population
   carries at least COVERAGE_THRESHOLD of the peak frequency. The permutation
   null is recomputed under the same definition.

2. The step 11 regression treats rsIDs as independent. They are not: 13
   PIK3R2 rsIDs come from one study and are in linkage disequilibrium, so
   p-values near 1e-16 are inflated by pseudo-replication. This step refits
   the models with cluster-robust (sandwich) standard errors, clustered by
   locus and, separately, by the study that contributes each variant.

Inputs:
  data_clean/final_table.csv
  data_clean/annotation_resolution_v3.csv
  variantAnnotations/study_parameters.tsv
  variantAnnotations/var_{drug,pheno,fa}_ann.tsv     (PMID per annotation)

Outputs (results/robustness/):
  coverage_by_stratum.csv        -- binary vs. continuous coverage per stratum
  coverage_per_rsid.csv          -- per-rsID coverage ratios
  coverage_permutation_null.csv  -- null under the coverage definition
  regression_clustered.csv       -- models with naive and clustered SEs
  coverage_clustering_report.txt

No new dependencies.
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

PRIMARY_POPS       = ["AFR", "AMR", "EAS", "NFE", "SAS"]
COVERAGE_THRESHOLD = 0.50     # studied pop must carry >=50% of peak AF
N_PERM             = 1_000
SEED               = 20261025

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


def norm_sf_two_sided(z):
    return math.erfc(abs(z) / math.sqrt(2.0))


# ── Logistic regression with optional cluster-robust SEs ────────────────────
def logistic_fit(X, y, names, clusters=None, ridge=1e-6, max_iter=200, tol=1e-11):
    """Fit logistic regression. If `clusters` is given, return cluster-robust
    (sandwich) standard errors alongside the naive ones."""
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    n, k = X.shape
    beta = np.zeros(k)
    for _ in range(max_iter):
        p = 1.0 / (1.0 + np.exp(-np.clip(X @ beta, -35, 35)))
        W = np.clip(p * (1 - p), 1e-10, None)
        H = (X * W[:, None]).T @ X + ridge * np.eye(k)
        grad = X.T @ (y - p) - ridge * beta
        try:
            step = np.linalg.solve(H, grad)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(H, grad, rcond=None)[0]
        beta_new = beta + step
        if np.max(np.abs(beta_new - beta)) < tol:
            beta = beta_new
            break
        beta = beta_new

    p = 1.0 / (1.0 + np.exp(-np.clip(X @ beta, -35, 35)))
    W = np.clip(p * (1 - p), 1e-10, None)
    bread = np.linalg.pinv((X * W[:, None]).T @ X + ridge * np.eye(k))
    se_naive = np.sqrt(np.maximum(np.diag(bread), 0))

    if clusters is None:
        se = se_naive
        n_clusters = np.nan
    else:
        clusters = np.asarray(clusters)
        u = y - p
        meat = np.zeros((k, k))
        uniq = pd.unique(clusters)
        for g in uniq:
            m = clusters == g
            sg = X[m].T @ u[m]
            meat += np.outer(sg, sg)
        G = len(uniq)
        # standard small-sample correction
        corr = (G / max(G - 1, 1)) * ((n - 1) / max(n - k, 1))
        V = bread @ (meat * corr) @ bread
        se = np.sqrt(np.maximum(np.diag(V), 0))
        n_clusters = G

    z = np.divide(beta, se, out=np.zeros_like(beta), where=se > 0)
    z_naive = np.divide(beta, se_naive, out=np.zeros_like(beta), where=se_naive > 0)
    return pd.DataFrame({
        "term": names,
        "coef": beta.round(4),
        "odds_ratio": np.exp(beta).round(3),
        "se_naive": se_naive.round(4),
        "p_naive": [norm_sf_two_sided(v) for v in z_naive],
        "se_used": se.round(4),
        "or_ci_low": np.exp(beta - 1.96 * se).round(3),
        "or_ci_high": np.exp(beta + 1.96 * se).round(3),
        "p_value": [norm_sf_two_sided(v) for v in z],
        "n_obs": n,
        "n_clusters": n_clusters,
    })


# ── Load and assemble the per-rsID table ────────────────────────────────────
say("=== CONTINUOUS COVERAGE METRIC AND CLUSTER-ROBUST INFERENCE ===")
say()

df = pd.read_csv(CLEAN / "final_table.csv", low_memory=False)
res = pd.read_csv(CLEAN / "annotation_resolution_v3.csv", low_memory=False)
sp = pd.read_csv(VARANN / "study_parameters.tsv", sep="\t", low_memory=False)

sp["aid"] = pd.to_numeric(sp["Variant Annotation ID"], errors="coerce")
res["aid"] = pd.to_numeric(res["annotation_id"], errors="coerce")
sp["group_label"] = (sp["Biogeographical Groups"].fillna("")
                     .astype(str).str.split(",").str[0].str.strip())
sp["mapped_pop"] = sp["group_label"].map(ANCESTRY_MAP)
sp["participants"] = (pd.to_numeric(sp["Study Cases"], errors="coerce").fillna(0)
                      + pd.to_numeric(sp["Study Controls"], errors="coerce").fillna(0))

studies = sp.merge(res[["aid", "rsID"]].drop_duplicates(), on="aid", how="inner")

# PMID per annotation, for study-level clustering
pmids = []
for stem in ["var_drug_ann", "var_pheno_ann", "var_fa_ann"]:
    t = pd.read_csv(VARANN / f"{stem}.tsv", sep="\t", low_memory=False)
    pmids.append(t[["Variant Annotation ID", "PMID"]])
pmids = pd.concat(pmids, ignore_index=True)
pmids["aid"] = pd.to_numeric(pmids["Variant Annotation ID"], errors="coerce")
pmid_by_rsid = (pmids.merge(res[["aid", "rsID"]].drop_duplicates(), on="aid", how="inner")
                .dropna(subset=["PMID"])
                .groupby("rsID")["PMID"]
                .agg(lambda s: s.value_counts().idxmax())
                .rename("primary_pmid"))

# Per-rsID frequency profile (first row per rsID; frequencies are variant-level)
freq_cols = [f"{p}_af" for p in PRIMARY_POPS]
per_rsid = (df[df["skew_assessable"]]
            .groupby("rsID")
            .agg({**{c: "first" for c in freq_cols},
                  "Gene": "first",
                  "max_AF": "first",
                  "max_pop": "first",
                  "skew_flag_5x": "any",
                  "pgkb_evidence_tier": lambda s: "high" if (s == "high").any() else "low"})
            .reset_index()
            .rename(columns={"pgkb_evidence_tier": "tier",
                             "skew_flag_5x": "high_skew"}))

studied = (studies.dropna(subset=["mapped_pop"])
           .groupby("rsID")["mapped_pop"].apply(lambda s: sorted(set(s)))
           .rename("studied_pops"))
n_studies = studies.groupby("rsID").size().rename("n_studies")
participants = studies.groupby("rsID")["participants"].sum().rename("total_participants")

t = (per_rsid.merge(studied, on="rsID", how="inner")
     .merge(n_studies, on="rsID", how="left")
     .merge(participants, on="rsID", how="left")
     .merge(pmid_by_rsid, on="rsID", how="left"))
t["locus"] = t["Gene"].astype(str).str.split(";").str[0].str.strip()
t["stratum"] = (np.where(t["high_skew"], "high skew", "low skew")
                + " / " + t["tier"] + " evidence")
say(f"rsIDs with frequency profile and >=1 mapped cohort: {len(t):,}")
say()


# ── 1. Continuous coverage metric ───────────────────────────────────────────
say(f"--- 1. Coverage metric (threshold = {COVERAGE_THRESHOLD:.0%} of peak AF) ---")


def coverage(row):
    pops = row["studied_pops"]
    peak = row["max_AF"]
    if not pops or pd.isna(peak) or peak <= 0:
        return pd.Series({"studied_max_af": np.nan, "coverage_ratio": np.nan})
    vals = [row[f"{p}_af"] for p in pops
            if f"{p}_af" in row.index and pd.notna(row[f"{p}_af"])]
    if not vals:
        return pd.Series({"studied_max_af": np.nan, "coverage_ratio": np.nan})
    smax = max(vals)
    return pd.Series({"studied_max_af": smax, "coverage_ratio": smax / peak})


t = pd.concat([t, t.apply(coverage, axis=1)], axis=1)
t["binary_excludes_peak"] = t.apply(
    lambda r: 0 if r["max_pop"] in r["studied_pops"] else 1, axis=1)
t["uncovered_50"] = (t["coverage_ratio"] < COVERAGE_THRESHOLD).astype(int)

cov = (t.dropna(subset=["coverage_ratio"])
       .groupby("stratum")
       .agg(rsids=("coverage_ratio", "size"),
            pct_binary_excludes_peak=("binary_excludes_peak",
                                      lambda s: round(100 * s.mean(), 1)),
            median_coverage_ratio=("coverage_ratio", lambda s: round(s.median(), 3)),
            pct_uncovered_50=("uncovered_50", lambda s: round(100 * s.mean(), 1)))
       .reset_index())
say(cov.to_string(index=False))
cov.to_csv(OUT / "coverage_by_stratum.csv", index=False)
t.drop(columns=["studied_pops"]).assign(
    studied_pops=t["studied_pops"].apply(";".join)
).to_csv(OUT / "coverage_per_rsid.csv", index=False)
say()
say("  The binary metric counts a variant as mismatched whenever the single")
say("  highest-frequency population was not studied, even if a studied")
say("  population carries a near-identical frequency. The coverage ratio")
say("  measures how much of the peak frequency the studied cohorts actually")
say("  reach, so near-ties no longer register as neglect.")
say()


# ── 2. Permutation null under the coverage definition ───────────────────────
say(f"--- 2. Permutation null under the coverage definition ({N_PERM:,} draws) ---")

glob_counts = t["studied_pops"].explode().value_counts()
pops = list(glob_counts.index)
weights = (glob_counts / glob_counts.sum()).reindex(pops).to_numpy()
weights = weights / weights.sum()
say("  rsID-level cohort marginals:")
for p, w in zip(pops, weights):
    say(f"    {p}: {100*w:.1f}%")

rng = np.random.default_rng(SEED)
sub = t.dropna(subset=["coverage_ratio"]).reset_index(drop=True)
af_mat = sub[[f"{p}_af" for p in pops]].to_numpy(float)
peak = sub["max_AF"].to_numpy(float)
ks = sub["studied_pops"].apply(len).to_numpy()

rows = []
for stratum, idx in sub.groupby("stratum").groups.items():
    idx = np.asarray(list(idx))
    obs_bin = 100 * sub.loc[idx, "binary_excludes_peak"].mean()
    obs_cov = 100 * sub.loc[idx, "uncovered_50"].mean()
    acc_bin = acc_cov = 0.0
    for _ in range(N_PERM):
        miss_b = miss_c = 0
        for i, k in zip(idx, ks[idx]):
            k_eff = min(int(k), len(pops))
            drawn = rng.choice(len(pops), size=k_eff, replace=False, p=weights)
            vals = af_mat[i, drawn]
            vals = vals[~np.isnan(vals)]
            pk = peak[i]
            if len(vals) == 0 or pk <= 0:
                miss_b += 1
                miss_c += 1
                continue
            if pops[int(np.nanargmax(af_mat[i]))] not in [pops[d] for d in drawn]:
                miss_b += 1
            if vals.max() / pk < COVERAGE_THRESHOLD:
                miss_c += 1
        acc_bin += miss_b / len(idx)
        acc_cov += miss_c / len(idx)
    rows.append({
        "stratum": stratum, "rsids": len(idx),
        "observed_pct_binary": round(obs_bin, 1),
        "null_pct_binary": round(100 * acc_bin / N_PERM, 1),
        "observed_pct_uncovered": round(obs_cov, 1),
        "null_pct_uncovered": round(100 * acc_cov / N_PERM, 1),
    })

perm = pd.DataFrame(rows)
perm["binary_obs_minus_null"] = (perm["observed_pct_binary"]
                                 - perm["null_pct_binary"]).round(1)
perm["coverage_obs_minus_null"] = (perm["observed_pct_uncovered"]
                                   - perm["null_pct_uncovered"]).round(1)
say()
say(perm.to_string(index=False))
perm.to_csv(OUT / "coverage_permutation_null.csv", index=False)
say()


# ── 3. Cluster-robust regression ────────────────────────────────────────────
say("--- 3. Logistic regression with cluster-robust standard errors ---")

reg = t.dropna(subset=["coverage_ratio", "n_studies"]).copy()
reg["log_cohorts"] = np.log1p(reg["n_studies"])
reg["log_participants"] = np.log1p(reg["total_participants"].fillna(0))
reg["high_skew_i"] = reg["high_skew"].astype(int)
reg["high_ev_i"] = (reg["tier"] == "high").astype(int)
reg["primary_pmid"] = reg["primary_pmid"].fillna("unknown").astype(str)

cols = ["high_skew_i", "high_ev_i", "log_cohorts", "log_participants"]
X = np.column_stack([np.ones(len(reg))] + [reg[c].to_numpy(float) for c in cols])
names = ["intercept"] + cols

fits = []
for outcome_name, yv in [("binary: peak population excluded",
                          reg["binary_excludes_peak"].to_numpy(float)),
                         (f"coverage: studied AF < {COVERAGE_THRESHOLD:.0%} of peak",
                          reg["uncovered_50"].to_numpy(float))]:
    for cl_name, cl in [("none (naive)", None),
                        ("locus (gene)", reg["locus"].to_numpy()),
                        ("study (primary PMID)", reg["primary_pmid"].to_numpy())]:
        f = logistic_fit(X, yv, names, clusters=cl)
        f.insert(0, "clustering", cl_name)
        f.insert(0, "outcome", outcome_name)
        fits.append(f)

    say(f"  Outcome: {outcome_name}   (n = {len(reg):,})")
    for cl_name in ["none (naive)", "locus (gene)", "study (primary PMID)"]:
        f = [x for x in fits if x["outcome"].iloc[0] == outcome_name
             and x["clustering"].iloc[0] == cl_name][0]
        ng = f["n_clusters"].iloc[0]
        tag = "" if pd.isna(ng) else f", {int(ng)} clusters"
        say(f"    clustering: {cl_name}{tag}")
        for _, r in f.iterrows():
            if r["term"] == "intercept":
                continue
            star = " *" if r["p_value"] < 0.05 else ""
            say(f"      {r['term']:<18} OR = {r['odds_ratio']:>7.3f}  "
                f"95% CI [{r['or_ci_low']:.3f}, {r['or_ci_high']:.3f}]  "
                f"p = {r['p_value']:.3g}{star}")
    say()

pd.concat(fits, ignore_index=True).to_csv(OUT / "regression_clustered.csv", index=False)

say("  Clustered SEs allow for shared frequency profiles and linkage within a")
say("  locus, and for several variants entering through the same publication.")
say("  Any term whose significance depends on the naive SEs was resting on")
say("  pseudo-replication.")
say()

report = OUT / "coverage_clustering_report.txt"
report.write_text("\n".join(lines), encoding="utf-8")
print(f"\nSaved report: {report}")
print(f"Saved tables to: {OUT}")
print("\nDone.")
