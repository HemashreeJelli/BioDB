from pathlib import Path
import re
import time
import requests
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent
RSID_FILE = BASE_DIR / "data_clean" / "rsid_list.txt"
MAP_FILE = BASE_DIR / "data_clean" / "rsid_to_variant_ids.csv"
CAND_FILE = BASE_DIR / "data_clean" / "gnomad_candidates.csv"

API_URL = "https://gnomad.broadinstitute.org/api"
DATASET = "gnomad_r4"
BATCH = 20
DELAY = 7
MAX_RETRIES = 4
POPS = ["afr", "amr", "eas", "nfe", "sas", "fin", "mid"]

RS_RE = re.compile(r"^rs\d+$")
VID_RE = re.compile(r"^[0-9XYM]+-\d+-[ACGT]+-[ACGT]+$")

session = requests.Session()


def post(query, variables):
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            r = session.post(API_URL, json={"query": query, "variables": variables}, timeout=90)
            if r.status_code == 429 or r.status_code >= 500:
                wait = DELAY * attempt * 3
                print(f"  HTTP {r.status_code}; waiting {wait}s (attempt {attempt})")
                time.sleep(wait)
                continue
            r.raise_for_status()
            return r.json()
        except (requests.RequestException, ValueError) as exc:
            wait = DELAY * attempt * 3
            print(f"  request problem: {exc}; waiting {wait}s (attempt {attempt})")
            time.sleep(wait)
    raise RuntimeError("API failed repeatedly. Progress is saved; re-run to resume.")


def failed_aliases(resp, n, prefix):
    """Aliases whose result cannot be trusted (GraphQL error on that alias, or whole response bad)."""
    if resp.get("data") is None:
        print("  whole batch failed:", str(resp.get("errors"))[:300])
        return {f"{prefix}{j}" for j in range(n)}
    bad = set()
    for e in resp.get("errors") or []:
        path = e.get("path")
        if path:
            bad.add(str(path[0]))
        else:
            print("  unattributed error:", str(e)[:300])
    return bad


def af(ac, an):
    return None if ac is None or an in (None, 0) else ac / an


def load_csv(path, columns):
    if path.exists():
        return pd.read_csv(path, dtype=str, keep_default_na=False)
    return pd.DataFrame(columns=columns)


def stage1_search(rsids):
    df = load_csv(MAP_FILE, ["rsID", "variant_id"])
    rows = df.to_dict("records")
    searched = set(df["rsID"])
    todo = [r for r in rsids if r not in searched]
    print(f"Stage 1: {len(todo)} rsIDs still to search")

    for i in range(0, len(todo), BATCH):
        chunk = todo[i:i + BATCH]
        parts = [f'a{j}: variant_search(query: "{rs}", dataset: $dataset) {{ variant_id }}'
                 for j, rs in enumerate(chunk)]
        resp = post("query($dataset: DatasetId!) {\n" + "\n".join(parts) + "\n}", {"dataset": DATASET})
        data = resp.get("data") or {}
        bad = failed_aliases(resp, len(chunk), "a")
        for j, rs in enumerate(chunk):
            alias = f"a{j}"
            result = data.get(alias)
            if alias in bad or result is None:
                continue  # not recorded -> retried on next run
            vids = [c["variant_id"] for c in result]
            if not vids:
                rows.append({"rsID": rs, "variant_id": ""})
            for v in vids:
                rows.append({"rsID": rs, "variant_id": v})
        pd.DataFrame(rows, columns=["rsID", "variant_id"]).to_csv(MAP_FILE, index=False)
        print(f"  searched {min(i + BATCH, len(todo))}/{len(todo)}")
        time.sleep(DELAY)


def stage2_fetch():
    mapping = load_csv(MAP_FILE, ["rsID", "variant_id"])
    wanted = sorted({v for v in mapping["variant_id"] if v})
    bad_format = [v for v in wanted if not VID_RE.match(v)]
    if bad_format:
        print("Skipping unexpected variant ID formats:", bad_format[:10])
    wanted = [v for v in wanted if VID_RE.match(v)]

    cols = (["variant_id", "fetch_status", "rsid_returned", "reference_genome", "chrom", "pos",
             "ref", "alt", "overall_ac", "overall_an", "overall_af", "dataset"]
            + [f"{p.upper()}_{m}" for p in POPS for m in ("ac", "an", "af")])
    done_df = load_csv(CAND_FILE, cols)
    rows = done_df.to_dict("records")
    done = set(done_df["variant_id"])
    todo = [v for v in wanted if v not in done]
    print(f"Stage 2: {len(todo)} candidate variants still to fetch")

    for i in range(0, len(todo), BATCH):
        chunk = todo[i:i + BATCH]
        parts = [f'''v{j}: variant(variantId: "{v}", dataset: $dataset) {{
          variant_id rsid reference_genome chrom pos ref alt
          joint {{ ac an populations {{ id ac an }} }} }}''' for j, v in enumerate(chunk)]
        resp = post("query($dataset: DatasetId!) {\n" + "\n".join(parts) + "\n}", {"dataset": DATASET})
        data = resp.get("data") or {}
        bad = failed_aliases(resp, len(chunk), "v")
        for j, vid in enumerate(chunk):
            alias = f"v{j}"
            if alias in bad:
                continue
            v = data.get(alias)
            if not v:
                rows.append({"variant_id": vid, "fetch_status": "no_data", "dataset": DATASET})
                continue
            joint = v.get("joint") or {}
            pops = {p["id"]: p for p in joint.get("populations", []) if p.get("id")}
            row = {"variant_id": vid, "fetch_status": "ok", "rsid_returned": v.get("rsid"),
                   "reference_genome": v.get("reference_genome"), "chrom": v.get("chrom"),
                   "pos": v.get("pos"), "ref": v.get("ref"), "alt": v.get("alt"),
                   "overall_ac": joint.get("ac"), "overall_an": joint.get("an"),
                   "overall_af": af(joint.get("ac"), joint.get("an")), "dataset": DATASET}
            for p in POPS:
                ac, an = pops.get(p, {}).get("ac"), pops.get(p, {}).get("an")
                row[f"{p.upper()}_ac"], row[f"{p.upper()}_an"] = ac, an
                row[f"{p.upper()}_af"] = af(ac, an)
            rows.append(row)
        pd.DataFrame(rows, columns=cols).to_csv(CAND_FILE, index=False)
        print(f"  fetched {min(i + BATCH, len(todo))}/{len(todo)}")
        time.sleep(DELAY)


def main():
    rsids = list(dict.fromkeys(
        l.strip() for l in RSID_FILE.read_text(encoding="utf-8").splitlines() if l.strip()))
    invalid = [r for r in rsids if not RS_RE.match(r)]
    if invalid:
        print("Dropping invalid rsIDs:", invalid[:10])
    rsids = [r for r in rsids if RS_RE.match(r)]
    print("Total rsIDs:", len(rsids))

    stage1_search(rsids)
    stage2_fetch()

    mapping = load_csv(MAP_FILE, ["rsID", "variant_id"])
    cand = load_csv(CAND_FILE, ["variant_id", "fetch_status"])
    searched = set(mapping["rsID"])
    print("\n=== SUMMARY ===")
    print("rsIDs searched:", len(searched), "of", len(rsids))
    per_rs = mapping[mapping["variant_id"] != ""].groupby("rsID")["variant_id"].count()
    print("rsIDs with >=1 candidate:", len(per_rs))
    print("rsIDs with NO candidate:", len(searched) - len(per_rs))
    print("Candidates per rsID:\n", per_rs.value_counts().sort_index().to_string())
    print("Fetch status:\n", cand["fetch_status"].value_counts().to_string())
    if len(searched) < len(rsids):
        print("\nSome rsIDs were not searched; just re-run the script to resume.")


if __name__ == "__main__":
    main()