from pathlib import Path
import re
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent
CLEAN = BASE_DIR / "data_clean"
VAR_DIR = BASE_DIR / "variantAnnotations"

RES_FILE = CLEAN / "annotation_resolution.csv"
CORE_FILE = CLEAN / "core_variants.csv"
OUT_RES = CLEAN / "annotation_resolution_v2.csv"
OUT_REVIEW = CLEAN / "manual_review_list_v2.csv"
OUT_CORE = CLEAN / "core_rows_usable_v2.csv"
OUT_DIAG = CLEAN / "core_rows_no_annotation_diagnostic.csv"

COMP = str.maketrans("ACGT", "TGCA")
SINGLE_STRAND = re.compile(r"single[- ]strand", re.I)
COMPLEMENTED = re.compile(r"complemented", re.I)
PLUS = re.compile(r"(plus|positive|forward|\+)\s*(chromosomal\s*)?strand", re.I)
RISK = re.compile(
    r"(minus|negative|reverse)\s*(chromosomal\s*)?strand|complement(ary)?\s+strand"
    r"|antisense|coding strand", re.I)
ANY = re.compile(r"strand|complement|antisense", re.I)
SNIPPET = re.compile(r"[^.;]*(strand|complement|antisense)[^.;]*", re.I)
USABLE = {"accepted", "accepted_plus_strand_assumed", "accepted_complemented"}


def find_col(df, name):
    return next((c for c in df.columns if c.lower().replace(" ", "") == name), None)


def strand_note_v2(text):
    t = SINGLE_STRAND.sub(" ", text)
    if not ANY.search(t):
        return ""
    if COMPLEMENTED.search(t):
        return "converted_to_plus"
    if PLUS.search(t):
        return "plus"
    if RISK.search(t):
        return "not_converted_or_unclear"
    return "other_mention"


def resolve_v2(row):
    st, sn = row["match_status"], row["strand_note_v2"]
    if st == "matched_direct":
        return "accepted"
    if st == "matched_strand_ambiguous":
        if sn in ("not_converted_or_unclear", "other_mention"):
            return "review_strand_note"
        return "accepted_plus_strand_assumed"
    if st == "matched_flipped":
        if sn in ("converted_to_plus", "plus"):
            return "review_flipped_contradicts_note"
        return "accepted_complemented"
    return row["resolution"]  # keep the earlier decision for everything else


def to_plus(letters, resolution):
    out = letters.translate(COMP) if resolution == "accepted_complemented" else letters
    return "".join(sorted(out))


def snippet(text):
    m = SNIPPET.search(text)
    return m.group(0).strip()[:140] if m else ""


def main():
    ann = pd.read_csv(RES_FILE, dtype=str, keep_default_na=False)
    core = pd.read_csv(CORE_FILE, dtype=str, keep_default_na=False)
    rs_col, gene_col = find_col(core, "rsid"), find_col(core, "gene")
    if rs_col is None:
        print("No rsID column in core file. Columns:", list(core.columns))
        return

    ann["strand_note_v2"] = ann["notes"].map(strand_note_v2)
    ann["resolution_v2"] = ann.apply(resolve_v2, axis=1)
    ann["effect_alleles_plus"] = ann.apply(
        lambda r: to_plus(r["effect_alleles"], r["resolution_v2"]), axis=1)
    ann["comparison_alleles_plus"] = ann.apply(
        lambda r: to_plus(r["comparison_alleles"], r["resolution_v2"]), axis=1)
    ann.to_csv(OUT_RES, index=False)

    print("=== resolution_v2 (in-scope annotations) ===")
    print(ann["resolution_v2"].value_counts().to_string())

    print("\n=== v1 -> v2 changes ===")
    changed = ann[ann["resolution"] != ann["resolution_v2"]]
    print(changed.groupby(["resolution", "resolution_v2"]).size().to_string()
          if len(changed) else "none")

    amb = ann[ann["match_status"] == "matched_strand_ambiguous"]
    print("\n=== strand_note_v2 among strand-ambiguous matches ===")
    print(amb["strand_note_v2"].replace("", "none").value_counts().to_string())

    risky = ann[ann["strand_note_v2"].isin(["not_converted_or_unclear", "other_mention"])]
    print(f"\n--- notes still needing a look: {len(risky)} annotations; distinct snippets ---")
    if len(risky):
        print(risky["notes"].map(snippet).value_counts().head(15).to_string())

    review = ann[ann["resolution_v2"].str.startswith("review")].copy()
    review["notes"] = review["notes"].str.slice(0, 300)
    review[["rsID", "gene", "drug", "alleles_raw", "comparison_raw", "site_alleles",
            "match_status", "resolution_v2", "strand_note_v2", "notes",
            "annotation_id", "source_file"]].to_csv(OUT_REVIEW, index=False)

    # core-row usability at (rsID, gene) level
    ok = ann["resolution_v2"].isin(USABLE)
    if gene_col:
        ann["_g"] = ann["gene"].str.strip()
        usable_pairs = set(zip(ann.loc[ok, "rsID"], ann.loc[ok, "_g"]))
        any_pairs = set(zip(ann["rsID"], ann["_g"]))

        def row_status(r):
            key = (r[rs_col], r[gene_col].strip())
            if key in usable_pairs:
                return "usable"
            if key in any_pairs:
                return "annotations_found_but_none_usable"
            return "no_matching_annotation_rows"
        core["usability"] = core.apply(row_status, axis=1)
    else:
        usable_rs = set(ann.loc[ok, "rsID"])
        core["usability"] = core[rs_col].map(lambda r: "usable" if r in usable_rs else "not_usable")
    core.to_csv(OUT_CORE, index=False)
    print("\n=== core rows (rsID x gene x drug) ===")
    print(core["usability"].value_counts().to_string())

    # diagnostic for core rows with no annotation rows at all
    if gene_col:
        frames = []
        for fname in ("var_drug_ann.tsv", "var_fa_ann.tsv", "var_pheno_ann.tsv"):
            path = VAR_DIR / fname
            if path.exists():
                d = pd.read_csv(path, sep="\t", dtype=str)
                if {"Variant/Haplotypes", "Gene"} <= set(d.columns):
                    frames.append(d[["Variant/Haplotypes", "Gene"]]
                                  .rename(columns={"Variant/Haplotypes": "rsID"}))
        genes_by_rs = {}
        if frames:
            allv = pd.concat(frames)
            genes_by_rs = (allv.groupby("rsID")["Gene"]
                           .agg(lambda s: ";".join(sorted(set(s.dropna())))).to_dict())
        missing = core[core["usability"] == "no_matching_annotation_rows"].copy()
        missing["rsid_in_var_files"] = missing[rs_col].map(lambda r: r in genes_by_rs)
        missing["genes_in_var_files"] = missing[rs_col].map(lambda r: genes_by_rs.get(r, ""))
        missing.to_csv(OUT_DIAG, index=False)
        print(f"\n=== core rows with no matching annotation rows: {len(missing)} ===")
        print("  rsID absent from all variant-annotation files:",
              int((~missing["rsid_in_var_files"]).sum()))
        print("  rsID present, but under a different gene string:",
              int(missing["rsid_in_var_files"].sum()))
        if len(missing):
            print(missing[[rs_col, gene_col, "genes_in_var_files"]]
                  .head(12).to_string(index=False))

    print("\nSaved:", OUT_RES, OUT_REVIEW, OUT_CORE, OUT_DIAG, sep="\n  ")


if __name__ == "__main__":
    main()