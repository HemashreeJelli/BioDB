"""
04_resolve_matches.py -- Step 2c: first-pass resolution of allele matches

Input : data_clean/annotation_allele_matches.csv, data_clean/core_variants.csv,
        variantAnnotations/var_*_ann.tsv (for the Notes column)
Output: data_clean/annotation_resolution.csv

Turns each annotation's match status into a decision (accepted / review_* /
excluded_*), using strand wording in the PharmGKB notes. Step 05 refines this.
"""

from pathlib import Path
import re
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent
CLEAN = BASE_DIR / "data_clean"
VAR_DIR = BASE_DIR / "variantAnnotations"

ANN_FILE = CLEAN / "annotation_allele_matches.csv"
CORE_FILE = CLEAN / "core_variants.csv"
OUT_RES = CLEAN / "annotation_resolution.csv"

COMP = str.maketrans("ACGT", "TGCA")
ANY_STRAND = re.compile(r"strand|complement|antisense", re.I)
PLUS = re.compile(r"(plus|positive|forward)( chromosomal)? strand", re.I)
MINUS = re.compile(r"(minus|negative|reverse)( chromosomal)? strand|antisense|coding strand", re.I)
SNIPPET = re.compile(r"[^.;]*(strand|complement|antisense)[^.;]*", re.I)
USABLE = {"accepted", "accepted_plus_strand_assumed"}


def find_col(df, name):
    return next((c for c in df.columns if c.lower().replace(" ", "") == name), None)


def strand_note_type(text):
    if not ANY_STRAND.search(text):
        return ""
    if MINUS.search(text):
        return "minus_or_mixed"
    if PLUS.search(text):
        return "plus"
    return "other_mention"


def parse_sites(s):
    return [set(x.strip()) for x in str(s).split("|") if x.strip()]


def resolve(row):
    st = row["match_status"]
    if st == "matched_direct":
        return "accepted"
    if st == "matched_strand_ambiguous":
        if row["strand_note"] in ("minus_or_mixed", "other_mention"):
            return "review_strand_note"
        return "accepted_plus_strand_assumed"
    if st == "matched_flipped":
        return "review_flipped"
    if st == "partial_allele_overlap":
        X = set(row["effect_alleles"]) | set(row["comparison_alleles"])
        for S in parse_sites(row["site_alleles"]):
            if all(a in S or a.translate(COMP) in S for a in X):
                return "review_mixed_strand_within_annotation"
        return "review_partial_overlap_other"
    if st == "allele_mismatch":
        return "review_mismatch"
    return "excluded_" + st


def main():
    ann = pd.read_csv(ANN_FILE, dtype=str, keep_default_na=False)
    core = pd.read_csv(CORE_FILE, dtype=str, keep_default_na=False)
    rs_col = find_col(core, "rsid")
    if rs_col is None:
        print("No rsID column in core file. Columns:", list(core.columns))
        return

    n_all = len(ann)
    ann = ann[ann["rsID"].isin(set(core[rs_col]))].copy()
    print(f"Annotation rows: {n_all} total, {len(ann)} in scope (rsID in core list)")

    parts = []
    for fname in ("var_drug_ann.tsv", "var_fa_ann.tsv", "var_pheno_ann.tsv"):
        path = VAR_DIR / fname
        if not path.exists():
            continue
        d = pd.read_csv(path, sep="\t", dtype=str)
        if {"Notes", "Variant Annotation ID"} <= set(d.columns):
            parts.append(pd.DataFrame({
                "source_file": fname,
                "annotation_id": d["Variant Annotation ID"],
                "notes": d["Notes"].fillna("")}))
    if parts:
        notes = pd.concat(parts).drop_duplicates(["source_file", "annotation_id"])
        ann = ann.merge(notes, how="left", on=["source_file", "annotation_id"])
    else:
        ann["notes"] = ""
    ann["notes"] = ann["notes"].fillna("")

    ann["strand_note"] = ann["notes"].map(strand_note_type)
    ann["allele_evidence"] = ann["n_alleles_total"].astype(int).map(
        lambda n: "two_or_more_alleles" if n >= 2 else "one_allele")
    ann["resolution"] = ann.apply(resolve, axis=1)
    ann.to_csv(OUT_RES, index=False)
    review = ann[ann["resolution"].str.startswith("review")]

    print("\n=== resolution (in-scope annotations) ===")
    print(ann["resolution"].value_counts().to_string())

    amb = ann[ann["match_status"] == "matched_strand_ambiguous"]
    print("\n=== strand notes among strand-ambiguous matches ===")
    print(amb["strand_note"].replace("", "none").value_counts().to_string())

    print("\n=== allele evidence among accepted annotations ===")
    print(ann[ann["resolution"].isin(USABLE)]["allele_evidence"].value_counts().to_string())

    odd = ann[ann["strand_note"].isin(["minus_or_mixed", "other_mention"])]
    print(f"\n--- strand notes needing a look: {len(odd)} annotations; distinct snippets ---")
    if len(odd):
        snips = odd["notes"].map(
            lambda t: (SNIPPET.search(t).group(0).strip() if SNIPPET.search(t) else "")[:160])
        print(snips.value_counts().head(12).to_string())

    for res in ("review_mixed_strand_within_annotation", "review_partial_overlap_other",
                "review_flipped", "review_mismatch"):
        sub = review[review["resolution"] == res]
        print(f"\n--- {res}: {len(sub)} annotations; first 8 ---")
        if len(sub):
            show = sub[["rsID", "gene", "alleles_raw", "comparison_raw", "site_alleles"]].head(8).astype(str)
            print(show.apply(lambda c: c.str.slice(0, 26)).to_string(index=False))

    # Core-row usability is decided later, in 06_core_usability.py.
    print("\nSaved:", OUT_RES)


if __name__ == "__main__":
    main()