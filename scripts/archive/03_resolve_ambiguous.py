
import json
import time
import requests

API_URL = "https://gnomad.broadinstitute.org/api"
DATASET = "gnomad_r4"

AMBIGUOUS_RSIDS = [
    "rs10012",
    "rs10123866",
    "rs10124893",
]

QUERY = """
query ResolveVariant($rsid: String!, $dataset: DatasetId!) {
  variant(rsid: $rsid, dataset: $dataset) {
    variant_id
    rsid
    reference_genome
    chrom
    pos
    ref
    alt
  }
}
"""

session = requests.Session()
session.headers.update({"Content-Type": "application/json"})

for rsid in AMBIGUOUS_RSIDS:
    print(f"\n{'=' * 60}")
    print(f"Checking {rsid}")

    response = session.post(
        API_URL,
        json={
            "query": QUERY,
            "variables": {
                "rsid": rsid,
                "dataset": DATASET,
            },
        },
        timeout=45,
    )

    print("HTTP status:", response.status_code)

    try:
        print(json.dumps(response.json(), indent=2))
    except ValueError:
        print(response.text[:2000])

    response.raise_for_status()
    time.sleep(7)
