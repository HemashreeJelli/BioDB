"""
05_refine_strand_notes.py -- Step 2d: refine resolution using strand notes

Input : data_clean/annotation_resolution.csv
Output: data_clean/annotation_resolution_v2.csv

Re-reads the PharmGKB notes with stricter strand rules, accepts flipped
matches as complemented when nothing contradicts that, and adds
effect_alleles_plus / comparison_alleles_plus (alleles on the + strand),
which step 08 uses to pick the annotated allele.
"""

from pathlib import Path
import re
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent
CLEAN = BASE_DIR / "data_clean"

RES_FILE = CLEAN / "annotation_resolution.csv"
OUT_RES = CLEAN / "annotation_resolution_v2.csv"

COMP = str.maketrans("ACGT", "TGCA")
SINGLE_STRAND = re.compile(r"single[- ]strand", re.I)
COMPLEMENTED = re.compile(r"complemented", re.I)
PLUS = re.compile(r"(plus|positive|forward|\+)\s*(chromosomal\s*)?strand", re.I)
RISK = re.compile(
    r"(minus|negative|reverse)\s*(chromosomal\s*)?strand|complement(ary)?\s+strand"
    r"|antisense|coding strand", re.I)
ANY = re.compile(r"strand|complement|antisense", re.I)
SNIPPET = re.compile(r"[^.;]*(strand|complement|antisense)[^.;]*", re.I)


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

    print("\nSaved:", OUT_RES)


if __name__ == "__main__":
    main()