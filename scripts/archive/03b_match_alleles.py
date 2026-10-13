from pathlib import Path
import re
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent
CLEAN = BASE_DIR / "data_clean"
VAR_DIR = BASE_DIR / "variantAnnotations"

CORE_FILE = CLEAN / "core_variants.csv"
MAP_FILE = CLEAN / "rsid_to_variant_ids.csv"
CAND_FILE = CLEAN / "gnomad_candidates.csv"
OUT_RS = CLEAN / "gnomad_matched_by_rsid.csv"
OUT_CORE = CLEAN / "core_with_frequencies.csv"
OUT_REVIEW = CLEAN / "gnomad_unresolved_review.csv"

POPS = ["afr", "amr", "eas", "nfe", "sas", "fin", "mid"]
COMP = str.maketrans("ACGT", "TGCA")


def comp(alleles):
    return {a.translate(COMP) for a in alleles}


def load_clinpgx_alleles():
    """rsID -> (set of single-letter alleles, set of non-standard tokens)."""
    comp_col = "Comparison Allele(s) or Genotype(s)"
    parts = []
    for fname in ("var_drug_ann.tsv", "var_fa_ann.tsv", "var_pheno_ann.tsv"):
        path = VAR_DIR / fname
        if not path.exists():
            print(f"WARNING: {fname} not found")
            continue
        d = pd.read_csv(path, sep="\t", dtype=str)
        if "Variant/Haplotypes" not in d.columns:
            print(f"WARNING: {fname} has no 'Variant/Haplotypes' column; columns = {list(d.columns)}")
            continue
        for col in ("Alleles", comp_col):
            if col in d.columns:
                parts.append(d[["Variant/Haplotypes", col]].set_axis(["rsID", "text"], axis=1))
            else:
                print(f"NOTE: {fname} has no column '{col}'")
    long = pd.concat(parts, ignore_index=True)
    long = long[long["rsID"].fillna("").str.match(r"^rs\d+$")]

    result = {}
    for rs, grp in long.groupby("rsID"):
        letters, nonstd = set(), set()
        for text in grp["text"].dropna():
            for tok in re.split(r"[^A-Za-z]+", str(text)):
                if not tok:
                    continue
                if re.fullmatch(r"[ACGT]{1,2}", tok.upper()):
                    letters.update(tok.upper())
                else:
                    nonstd.add(tok)
        result[rs] = (letters, nonstd)
    return result


def classify(S, cands):
    """Return (status, matched_candidate_or_None)."""
    if not cands:
        return "no_gnomad_candidate", None
    if not S:
        return "no_parsable_clinpgx_alleles", None
    snvs = [c for c in cands if len(c["ref"]) == 1 and len(c["alt"]) == 1]
    if not snvs:
        return "non_snv_unassessed", None
    if len(S) > 2:
        return "clinpgx_alleles_gt2", None

    cS = comp(S)
    if len(S) == 1:
        direct = [c for c in snvs if S <= {c["ref"], c["alt"]}]
        hits = direct or [c for c in snvs if cS <= {c["ref"], c["alt"]}]
        if len(hits) == 1:
            return "matched_single_allele_unverified", hits[0]
        return ("multiple_matches" if hits else "allele_mismatch"), None

    hits = []
    for c in snvs:
        a = {c["ref"], c["alt"]}
        d, x = S <= a, cS <= a
        if d and x:
            hits.append((c, "matched_strand_ambiguous"))
        elif d:
            hits.append((c, "matched_direct"))
        elif x:
            hits.append((c, "matched_strand_flipped"))
    if not hits:
        return "allele_mismatch", None
    if len(hits) > 1:
        return "multiple_matches", None
    return hits[0][1], hits[0][0]


def ref_af_approx(matched, cands):
    """1 - sum of ALT AFs of fetched candidates at the same site (approximation:
    AN differs slightly between variants, and other alts may not be listed)."""
    same = [c for c in cands if c["chrom"] == matched["chrom"] and c["pos"] == matched["pos"]]
    out = {}
    for p in POPS:
        afs = [c[f"{p.upper()}_af"] for c in same]
        out[f"{p.upper()}_ref_af_approx"] = (
            None if any(pd.isna(a) for a in afs) else max(0.0, 1 - sum(afs))
        )
    return out


def main():
    alleles = load_clinpgx_alleles()
    mapping = pd.read_csv(MAP_FILE, dtype=str, keep_default_na=False)
    cand = pd.read_csv(CAND_FILE, dtype=str, keep_default_na=False)
    cand = cand[cand["fetch_status"] == "ok"].copy()

    numeric = ["pos", "overall_ac", "overall_an", "overall_af"] + [
        f"{p.upper()}_{m}" for p in POPS for m in ("ac", "an", "af")
    ]
    for col in numeric:
        cand[col] = pd.to_numeric(cand[col], errors="coerce")
    cand_by_vid = cand.drop_duplicates("variant_id").set_index("variant_id").to_dict("index")

    vids_by_rs = {
        rs: [v for v in grp["variant_id"] if v]
        for rs, grp in mapping.groupby("rsID", sort=False)
    }

    rows, review = [], []
    for rs, vids in vids_by_rs.items():
        cands = [dict(cand_by_vid[v], variant_id=v) for v in vids if v in cand_by_vid]
        S, nonstd = alleles.get(rs, (set(), set()))
        if vids and not cands:
            status, matched = "candidates_not_fetched", None
        else:
            status, matched = classify(S, cands)

        row = {
            "rsID": rs,
            "match_status": status,
            "n_candidates": len(cands),
            "n_positions": len({(c["chrom"], c["pos"]) for c in cands}),
            "clinpgx_alleles": "".join(sorted(S)),
            "clinpgx_nonstandard_tokens": ";".join(sorted(nonstd)),
            "all_candidates": "; ".join(
                f"{c['variant_id']}[AC={int(c['overall_ac']) if pd.notna(c['overall_ac']) else 'NA'}]"
                for c in cands
            ),
        }
        if matched is not None:
            row.update({
                "variant_id": matched["variant_id"],
                "reference_genome": matched["reference_genome"],
                "chrom": matched["chrom"], "pos": matched["pos"],
                "ref": matched["ref"], "alt": matched["alt"],
                "overall_ac": matched["overall_ac"], "overall_an": matched["overall_an"],
                "overall_af": matched["overall_af"],
            })
            for p in POPS:
                for m in ("ac", "an", "af"):
                    row[f"{p.upper()}_{m}"] = matched[f"{p.upper()}_{m}"]
            row.update(ref_af_approx(matched, cands))
        rows.append(row)
        if not status.startswith("matched"):
            review.append({k: row[k] for k in (
                "rsID", "match_status", "clinpgx_alleles", "clinpgx_nonstandard_tokens",
                "n_candidates", "all_candidates")})

    rs_table = pd.DataFrame(rows)
    rs_table.to_csv(OUT_RS, index=False)
    pd.DataFrame(review).to_csv(OUT_REVIEW, index=False)

    # left join onto the ClinPGx core list (keeps every core row)
    core = pd.read_csv(CORE_FILE, dtype=str, keep_default_na=False)
    rs_col = next((c for c in core.columns if c.lower().replace(" ", "") == "rsid"), None)
    if rs_col is None:
        print("Could not find rsID column in core_variants.csv. Columns:", list(core.columns))
        return
    merged = core.merge(rs_table.rename(columns={"rsID": "_rs"}),
                        how="left", left_on=rs_col, right_on="_rs").drop(columns="_rs")
    merged["match_status"] = merged["match_status"].fillna("not_searched")
    merged.to_csv(OUT_CORE, index=False)

    print("=== rsID-level match status ===")
    print(rs_table["match_status"].value_counts().to_string())
    print("\n=== core-list rows (rsID x gene x drug) by status ===")
    print(merged["match_status"].value_counts().to_string())

    for status in ("allele_mismatch", "multiple_matches", "clinpgx_alleles_gt2"):
        sub = rs_table[rs_table["match_status"] == status].head(8)
        if len(sub):
            print(f"\n--- examples: {status} ---")
            print(sub[["rsID", "clinpgx_alleles", "all_candidates"]].to_string(index=False))

    print("\n--- rsIDs with >=4 candidates ---")
    big = rs_table[rs_table["n_candidates"] >= 4]
    print(big[["rsID", "n_candidates", "n_positions", "match_status"]].to_string(index=False))

    print("\n--- sanity check: the three investigated rsIDs ---")
    chk = rs_table[rs_table["rsID"].isin(["rs10012", "rs10123866", "rs10124893"])]
    print(chk[["rsID", "match_status", "clinpgx_alleles", "variant_id", "ref", "alt",
               "AFR_af", "EAS_af", "NFE_af"]].to_string(index=False))

    print("\nNon-standard allele tokens present for",
          int((rs_table["clinpgx_nonstandard_tokens"] != "").sum()), "rsIDs")
    print("\nSaved:", OUT_RS, OUT_CORE, OUT_REVIEW, sep="\n  ")


if __name__ == "__main__":
    main()