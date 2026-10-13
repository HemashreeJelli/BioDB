
import requests
import json
import time

API_URL = "https://gnomad.broadinstitute.org/api"

QUERY = """
query TestVariant($rsid: String!, $dataset: DatasetId!) {
  variant(rsid: $rsid, dataset: $dataset) {
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
      populations {
        id
        ac
        an
      }
    }
  }
}
"""

test_rsids = ["rs4149056", "rs9923231", "rs1801133"]

for rsid in test_rsids:
    print(f"\nTesting {rsid}...")

    try:
        response = requests.post(
            API_URL,
            json={
                "query": QUERY,
                "variables": {
                    "rsid": rsid,
                    "dataset": "gnomad_r4",
                },
            },
            headers={"Content-Type": "application/json"},
            timeout=30,
        )

        print("HTTP status:", response.status_code)
        response.raise_for_status()

        result = response.json()

        if result.get("errors"):
            print("GraphQL errors:")
            print(json.dumps(result["errors"], indent=2))
        else:
            variant = result.get("data", {}).get("variant")

            if variant is None:
                print("No variant returned.")
            else:
                print(json.dumps(variant, indent=2))

    except requests.RequestException as exc:
        print("Request failed:", exc)

    time.sleep(7)
