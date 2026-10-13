from pathlib import Path
import re
import time
import requests
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent
ERR_FILE = BASE_DIR / "data_clean" / "gnomad_errors.csv"
OUT_FILE = BASE_DIR / "data_clean" / "candidates_pilot.csv"

API_URL = "https://gnomad.broadinstitute.org/api"
DATASET = "gnomad_r4"
PILOT_N = 20
BATCH = 20          # forum reports a cap of 25 aliases per request
DELAY = 7
POPS = ["afr", "amr", "eas", "nfe", "sas", "fin", "mid"]

RS_RE = re.compile(r"^rs\d+$")
VID_RE = re.compile(r"^[0-9XYM]+-\d+-[ACGT]+-[ACGT]+$")

def post(query, variables):
    r = requests.post(API_URL, json={"query": query, "variables": variables}, timeout=60)
    print("HTTP", r.status_code)
    r.raise_for_status()
    return r.json()

def search_batch(rsids):
    parts = [f'a{i}: variant_search(query: "{rs}", dataset: $dataset) {{ variant_id }}'
             for i, rs in enumerate(rsids)]
    q = "query($dataset: DatasetId!) {\n" + "\n".join(parts) + "\n}"
    resp = post(q, {"dataset": DATASET})
    data = resp.get("data") or {}
    out = {}
    for i, rs in enumerate(rsids):
        out[rs] = [c["variant_id"] for c in (data.get(f"a{i}") or [])]
    if resp.get("errors"):
        print("SEARCH ERRORS:", resp["errors"])
    return out

def fetch_batch(vids):
    parts = [f'''v{i}: variant(variantId: "{v}", dataset: $dataset) {{
      variant_id rsid reference_genome chrom pos ref alt
      joint {{ ac an populations {{ id ac an }} }} }}'''
             for i, v in enumerate(vids)]
    q = "query($dataset: DatasetId!) {\n" + "\n".join(parts) + "\n}"
    resp = post(q, {"dataset": DATASET})
    data = resp.get("data") or {}
    if resp.get("errors"):
        print("FETCH ERRORS:", resp["errors"])
    return [data.get(f"v{i}") for i in range(len(vids))]

def af(ac, an):
    return None if ac is None or an in (None, 0) else ac / an

errs = pd.read_csv(ERR_FILE, dtype=str, keep_default_na=False)
rsids = errs.loc[errs["status"] == "ambiguous_variant", "rsID"].tolist()[:PILOT_N]
rsids = [r for r in rsids if RS_RE.match(r)]
print("Pilot rsIDs:", len(rsids))

cand_map = search_batch(rsids)
time.sleep(DELAY)

vids = sorted({v for lst in cand_map.values() for v in lst if VID_RE.match(v)})
print("Unique candidate variant IDs:", len(vids))

records = []
for start in range(0, len(vids), BATCH):
    chunk = vids[start:start + BATCH]
    variants = fetch_batch(chunk)
    time.sleep(DELAY)
    for vid, v in zip(chunk, variants):
        if not v:
            records.append({"variant_id": vid, "fetch_status": "no_data"})
            continue
        joint = v.get("joint") or {}
        pops = {p["id"]: p for p in joint.get("populations", [])}
        row = {"variant_id": vid, "fetch_status": "ok", "rsid_returned": v.get("rsid"),
               "chrom": v["chrom"], "pos": v["pos"], "ref": v["ref"], "alt": v["alt"],
               "overall_ac": joint.get("ac"), "overall_an": joint.get("an")}
        for p in POPS:
            ac, an = pops.get(p, {}).get("ac"), pops.get(p, {}).get("an")
            row[f"{p.upper()}_ac"], row[f"{p.upper()}_an"] = ac, an
            row[f"{p.upper()}_af"] = af(ac, an)
        records.append(row)

by_vid = {r["variant_id"]: r for r in records}
rows = []
for rs, lst in cand_map.items():
    if not lst:
        rows.append({"rsID": rs, "variant_id": None, "fetch_status": "no_candidates"})
    for vid in lst:
        rows.append({"rsID": rs, **by_vid.get(vid, {"variant_id": vid, "fetch_status": "not_fetched"})})

df = pd.DataFrame(rows)
df.to_csv(OUT_FILE, index=False)

print("\nCandidates per rsID:")
print(df.groupby("rsID")["variant_id"].count().value_counts().sort_index().to_string())
print("\nFetch status:")
print(df["fetch_status"].value_counts().to_string())
print("\nSaved:", OUT_FILE)
print(df[["rsID", "variant_id", "ref", "alt", "overall_ac", "overall_an"]].head(30).to_string(index=False))