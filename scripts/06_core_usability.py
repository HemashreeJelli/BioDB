"""
06_core_usability.py -- Step 2e: final resolution and core-row usability

Input : data_clean/annotation_resolution_v2.csv, data_clean/core_variants.csv
Output: data_clean/annotation_resolution_v3.csv
        data_clean/core_rows_usable_v3.csv     (core rows + usability column)
        data_clean/manual_review_list_v3.csv   (annotations still marked review_*)

Clears "(+) chromosomal strand" notes that step 05 flagged by mistake, then
marks a core row usable when any accepted annotation for that rsID shares a
gene with it.
"""

from pathlib import Path
import re
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent
CLEAN = BASE_DIR / "data_clean"
RES_IN = CLEAN / "annotation_resolution_v2.csv"
CORE_IN = CLEAN / "core_variants.csv"
RES_OUT = CLEAN / "annotation_resolution_v3.csv"
CORE_OUT = CLEAN / "core_rows_usable_v3.csv"
REVIEW_OUT = CLEAN / "manual_review_list_v3.csv"

USABLE = {"accepted", "accepted_plus_strand_assumed", "accepted_complemented"}
PLUS_PAREN = re.compile(r"\(\+\)\s*(?:chromosomal\s*)?strand", re.I)
OTHER_STRAND = re.compile(r"minus|reverse|negative|complement|antisense", re.I)


def find_col(df, name):
    return next((c for c in df.columns if c.lower().replace(" ", "") == name), None)


def genes(s):
    return {g.strip().upper() for g in re.split(r"[;,]", str(s)) if g.strip()}


ann = pd.read_csv(RES_IN, dtype=str, keep_default_na=False)
core = pd.read_csv(CORE_IN, dtype=str, keep_default_na=False)
rs_col, gene_col = find_col(core, "rsid"), find_col(core, "gene")

# Fix "(+) chromosomal strand" false alarms
fix = (ann["resolution_v2"] == "review_strand_note") \
    & ann["notes"].str.contains(PLUS_PAREN) \
    & ~ann["notes"].str.contains(OTHER_STRAND)
ann["resolution_v3"] = ann["resolution_v2"]
ann.loc[fix, "resolution_v3"] = "accepted_plus_strand_assumed"
ann.loc[fix, "strand_note_v2"] = "plus"
print("Reclassified '(+) strand' notes:", int(fix.sum()))
ann.to_csv(RES_OUT, index=False)

print("\n=== resolution_v3 ===")
print(ann["resolution_v3"].value_counts().to_string())

# Core usability with gene-token overlap
ann["_gset"] = ann["gene"].map(genes)
ann["_ok"] = ann["resolution_v3"].isin(USABLE)
by_rs = {rs: list(zip(g["_gset"], g["_ok"])) for rs, g in ann.groupby("rsID")}


def status(r):
    cands = [ok for gs, ok in by_rs.get(r[rs_col], []) if gs & genes(r[gene_col])]
    if not cands:
        return "no_matching_annotation_rows"
    return "usable" if any(cands) else "annotations_found_but_none_usable"


core["usability"] = core.apply(status, axis=1)
core.to_csv(CORE_OUT, index=False)
print("\n=== core rows (token-overlap gene matching) ===")
print(core["usability"].value_counts().to_string())
print(core.loc[core["usability"] == "no_matching_annotation_rows",
               [rs_col, gene_col]].head(10).to_string(index=False))

review = ann[ann["resolution_v3"].str.startswith("review")].copy()
review["notes"] = review["notes"].str.slice(0, 300)
review.sort_values(["resolution_v3", "rsID"]).to_csv(REVIEW_OUT, index=False)
print("\nReview rows:", len(review), "->", REVIEW_OUT)