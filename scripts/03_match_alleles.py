"""
03_match_alleles.py -- Step 2b: match annotation alleles to gnomAD sites

Input : data_clean/gnomad_candidates.csv, data_clean/rsid_to_variant_ids.csv,
        data_clean/core_variants.csv, variantAnnotations/var_*_ann.tsv
Output: data_clean/site_allele_frequencies.csv    (REF + ALT rows per site)
        data_clean/annotation_allele_matches.csv  (one row per variant annotation)
        data_clean/rsid_allele_match_summary.csv  (roll-up per rsID)
        data_clean/core_with_match_rollup.csv     (core rows + rsID roll-up)

Each annotation's alleles are checked against the gnomAD alleles directly and
as their complement (other strand). Only SNVs are assessed.
"""

from pathlib import Path
import re
import numpy as np
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent
CLEAN = BASE_DIR / "data_clean"
VAR_DIR = BASE_DIR / "variantAnnotations"

MAP_FILE = CLEAN / "rsid_to_variant_ids.csv"
CAND_FILE = CLEAN / "gnomad_candidates.csv"
CORE_FILE = CLEAN / "core_variants.csv"

OUT_ALLELES = CLEAN / "site_allele_frequencies.csv"
OUT_ANN = CLEAN / "annotation_allele_matches.csv"
OUT_RS = CLEAN / "rsid_allele_match_summary.csv"
OUT_CORE = CLEAN / "core_with_match_rollup.csv"

POPS = ["afr", "amr", "eas", "nfe", "sas", "fin", "mid"]
GROUPS = ["overall"] + POPS
METRICS = ("ac", "an", "af")
COMP = str.maketrans("ACGT", "TGCA")
COMP_COL = "Comparison Allele(s) or Genotype(s)"
MATCHED = {"matched_direct", "matched_flipped", "matched_strand_ambiguous"}


def col(group, metric):
    return f"overall_{metric}" if group == "overall" else f"{group.upper()}_{metric}"


def comp(alleles):
    return {a.translate(COMP) for a in alleles}


def parse_alleles(text):
    letters, nonstd = set(), set()
    if text is None or (isinstance(text, float) and np.isnan(text)):
        return letters, nonstd
    for tok in re.split(r"[^A-Za-z]+", str(text)):
        if not tok:
            continue
        if re.fullmatch(r"[ACGT]{1,2}", tok.upper()):
            letters.update(tok.upper())
        else:
            nonstd.add(tok)
    return letters, nonstd


def match(X, sites):
    cX = comp(X)
    hits = []
    for s in sites:
        direct = X <= s["alleles"]
        flipped = cX <= s["alleles"]
        if direct and flipped:
            hits.append((s, "matched_strand_ambiguous"))
        elif direct:
            hits.append((s, "matched_direct"))
        elif flipped:
            hits.append((s, "matched_flipped"))
    if len(hits) == 1:
        return hits[0][1], hits[0][0]
    if len(hits) > 1:
        return "multiple_sites", None
    overlap = any((X & s["alleles"]) or (cX & s["alleles"]) for s in sites)
    return ("partial_allele_overlap" if overlap else "allele_mismatch"), None


def main():
    # ---------- candidates ----------
    cand = pd.read_csv(CAND_FILE, dtype=str, keep_default_na=False)
    cand = cand[cand["fetch_status"] == "ok"].copy()
    for c in ["pos"] + [col(g, m) for g in GROUPS for m in METRICS]:
        cand[c] = pd.to_numeric(cand[c], errors="coerce")
    cand = cand.drop_duplicates("variant_id")
    cand_by_vid = {r["variant_id"]: r for r in cand.to_dict("records")}

    mapping = pd.read_csv(MAP_FILE, dtype=str, keep_default_na=False)

    # ---------- site-level allele table ----------
    sites_by_rs, status_pre, allele_rows = {}, {}, []
    for rs, grp in mapping.groupby("rsID", sort=False):
        vids = [v for v in grp["variant_id"] if v]
        if not vids:
            status_pre[rs] = "no_gnomad_candidate"
            continue
        rows = [cand_by_vid[v] for v in vids if v in cand_by_vid]
        if not rows:
            status_pre[rs] = "candidates_not_fetched"
            continue
        snvs = [r for r in rows if len(r["ref"]) == 1 and len(r["alt"]) == 1]
        if not snvs:
            status_pre[rs] = "non_snv_unassessed"
            continue

        by_site = {}
        for r in snvs:
            by_site.setdefault((r["chrom"], int(r["pos"])), []).append(r)

        sites = []
        for (chrom, pos), srows in by_site.items():
            if len({r["ref"] for r in srows}) != 1:
                print("WARNING: inconsistent REF at", rs, chrom, pos)
            ref = srows[0]["ref"]
            sites.append({"chrom": chrom, "pos": pos,
                          "alleles": {ref} | {r["alt"] for r in srows}})

            ref_row = {"rsID": rs, "chrom": chrom, "pos": pos, "allele": ref,
                       "allele_type": "REF", "variant_id": "",
                       "af_method": "1 - sum(ALT AF), approximate",
                       "n_alt_records": len(srows)}
            for g in GROUPS:
                afs = [r[col(g, "af")] for r in srows]
                valid_an = [r[col(g, "an")] for r in srows if pd.notna(r[col(g, "an")])]
                ref_row[col(g, "ac")] = np.nan
                ref_row[col(g, "an")] = min(valid_an) if valid_an else np.nan
                ref_row[col(g, "af")] = (np.nan if any(pd.isna(a) for a in afs)
                                         else max(0.0, 1 - sum(afs)))
            allele_rows.append(ref_row)

            for r in srows:
                row = {"rsID": rs, "chrom": chrom, "pos": pos, "allele": r["alt"],
                       "allele_type": "ALT", "variant_id": r["variant_id"],
                       "af_method": "AC/AN", "n_alt_records": len(srows)}
                for g in GROUPS:
                    for m in METRICS:
                        row[col(g, m)] = r[col(g, m)]
                allele_rows.append(row)
        sites_by_rs[rs] = sites

    pd.DataFrame(allele_rows).to_csv(OUT_ALLELES, index=False)

    # ---------- ClinPGx annotations (one row each, NOT pooled) ----------
    parts = []
    for fname in ("var_drug_ann.tsv", "var_fa_ann.tsv", "var_pheno_ann.tsv"):
        path = VAR_DIR / fname
        if not path.exists():
            print(f"WARNING: {fname} not found")
            continue
        d = pd.read_csv(path, sep="\t", dtype=str)
        if "Variant/Haplotypes" not in d.columns:
            print(f"WARNING: {fname} lacks 'Variant/Haplotypes'; columns = {list(d.columns)}")
            continue
        d = d[d["Variant/Haplotypes"].fillna("").str.match(r"^rs\d+$")].copy()
        parts.append(pd.DataFrame({
            "source_file": fname,
            "annotation_id": d.get("Variant Annotation ID"),
            "rsID": d["Variant/Haplotypes"],
            "gene": d.get("Gene"),
            "drug": d.get("Drug(s)"),
            "alleles_raw": d.get("Alleles"),
            "comparison_raw": d.get(COMP_COL),
        }))
    ann = pd.concat(parts, ignore_index=True)
    print("Annotation rows with rsIDs:", len(ann))

    results = []
    for r in ann.to_dict("records"):
        E, nE = parse_alleles(r["alleles_raw"])
        K, nK = parse_alleles(r["comparison_raw"])
        X = E | K
        rs = r["rsID"]
        sites = sites_by_rs.get(rs)
        chrom = pos = None
        if rs in status_pre:
            status = status_pre[rs]
        elif sites is None:
            status = "rsid_not_in_candidate_map"
        elif not X:
            status = "no_parsable_alleles"
        else:
            status, site = match(X, sites)
            if site is not None:
                chrom, pos = site["chrom"], site["pos"]
        results.append({
            **r,
            "effect_alleles": "".join(sorted(E)),
            "comparison_alleles": "".join(sorted(K)),
            "n_alleles_total": len(X),
            "nonstandard_tokens": ";".join(sorted(nE | nK)),
            "match_status": status,
            "matched_chrom": chrom, "matched_pos": pos,
            "site_alleles": " | ".join("".join(sorted(s["alleles"])) for s in sites) if sites else "",
        })
    ann_df = pd.DataFrame(results)
    ann_df.to_csv(OUT_ANN, index=False)

    # ---------- rsID roll-up ----------
    counts = ann_df.groupby(["rsID", "match_status"]).size().unstack(fill_value=0)
    n_total = counts.sum(axis=1)
    mcols = [c for c in counts.columns if c in MATCHED]
    n_matched = counts[mcols].sum(axis=1) if mcols else n_total * 0
    rs_summary = counts.copy()
    rs_summary["n_annotations"] = n_total
    rs_summary["n_matched"] = n_matched
    rs_summary["rsid_rollup"] = np.where(
        n_matched == n_total, "all_annotations_matched",
        np.where(n_matched > 0, "partly_matched", "none_matched"))
    rs_summary = rs_summary.reset_index()
    rs_summary.to_csv(OUT_RS, index=False)

    # ---------- report ----------
    print("\n=== annotation-level match status ===")
    print(ann_df["match_status"].value_counts().to_string())
    print("\n=== rsID-level roll-up ===")
    print(rs_summary["rsid_rollup"].value_counts().to_string())

    matched = ann_df[ann_df["match_status"].isin(MATCHED)]
    print("\nMatched annotations that list only ONE allele letter (weaker evidence):",
          int((matched["n_alleles_total"] == 1).sum()), "of", len(matched))

    def trunc(df):
        df = df.copy()
        for c in df.columns:
            df[c] = df[c].astype(str).str.slice(0, 28)
        return df

    cols = ["rsID", "gene", "drug", "alleles_raw", "comparison_raw", "site_alleles"]
    for status in ("partial_allele_overlap", "allele_mismatch", "multiple_sites"):
        sub = ann_df[ann_df["match_status"] == status]
        print(f"\n--- examples: {status} ({len(sub)} annotations) ---")
        if len(sub):
            print(trunc(sub[cols].head(10)).to_string(index=False))

    # ---------- core-row level ----------
    core = pd.read_csv(CORE_FILE, dtype=str, keep_default_na=False)
    rs_col = next((c for c in core.columns if c.lower().replace(" ", "") == "rsid"), None)
    if rs_col:
        merged = core.merge(
            rs_summary[["rsID", "rsid_rollup"]].rename(columns={"rsID": "_rs"}),
            how="left", left_on=rs_col, right_on="_rs").drop(columns="_rs")
        merged["rsid_rollup"] = merged["rsid_rollup"].fillna("no_variant_annotation_rows")
        merged.to_csv(OUT_CORE, index=False)
        print("\n=== core rows by roll-up ===")
        print(merged["rsid_rollup"].value_counts().to_string())

    print("\nSaved:", OUT_ALLELES, OUT_ANN, OUT_RS, OUT_CORE, sep="\n  ")


if __name__ == "__main__":
    main()