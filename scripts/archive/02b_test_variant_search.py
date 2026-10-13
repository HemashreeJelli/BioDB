from pathlib import Path
import json
import time
import requests

BASE_DIR = Path(__file__).resolve().parent.parent
OUT = BASE_DIR / "data_clean" / "test_candidates.json"

API_URL = "https://gnomad.broadinstitute.org/api"
DATASET = "gnomad_r4"
TEST_RSID = "rs10012"      # CYP1B1, ClinPGx alleles: C vs G
DELAY = 7

SEARCH_QUERY = """
query Search($q: String!, $dataset: DatasetId!) {
  variant_search(query: $q, dataset: $dataset) {
    variant_id
  }
}
"""

VARIANT_QUERY = """
query Var($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    variant_id
    rsid
    reference_genome
    chrom
    pos
    ref
    alt
    joint {
      ac
      an
      populations { id ac an }
    }
  }
}
"""

def post(query, variables):
    r = requests.post(API_URL, json={"query": query, "variables": variables}, timeout=45)
    print("HTTP", r.status_code)
    return r.json()

search = post(SEARCH_QUERY, {"q": TEST_RSID, "dataset": DATASET})
print("\nRAW SEARCH RESPONSE:")
print(json.dumps(search, indent=2))

candidates = (search.get("data") or {}).get("variant_search") or []
ids = [c["variant_id"] for c in candidates]
print("\nCandidate variant IDs:", ids)

results = {"rsid": TEST_RSID, "search": search, "candidates": []}

for vid in ids:
    time.sleep(DELAY)
    resp = post(VARIANT_QUERY, {"variantId": vid, "dataset": DATASET})
    v = (resp.get("data") or {}).get("variant")
    results["candidates"].append(resp)
    if not v:
        print(vid, "-> no data / error:", resp.get("errors"))
        continue
    joint = v.get("joint") or {}
    print(f"\n{vid}: rsid={v.get('rsid')} {v['chrom']}:{v['pos']} "
          f"REF={v['ref']} ALT={v['alt']} AC={joint.get('ac')} AN={joint.get('an')}")
    for p in joint.get("populations", []):
        if p["id"] in ("afr", "amr", "eas", "nfe", "sas", "fin", "mid"):
            print(f"   {p['id']}: AC={p['ac']} AN={p['an']}")

OUT.write_text(json.dumps(results, indent=2), encoding="utf-8")
print("\nSaved to", OUT)