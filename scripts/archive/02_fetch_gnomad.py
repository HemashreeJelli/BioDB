
from pathlib import Path
import json
import time
import requests
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent
RSID_FILE = BASE_DIR / "data_clean" / "rsid_list.txt"
OUTPUT_DIR = BASE_DIR / "data_clean"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_FILE = OUTPUT_DIR / "gnomad_frequencies.csv"
ERROR_FILE = OUTPUT_DIR / "gnomad_errors.csv"

API_URL = "https://gnomad.broadinstitute.org/api"
DATASET = "gnomad_r4"

REQUEST_DELAY = 7
MAX_RETRIES = 3
SAVE_EVERY = 25

POPULATIONS = ["afr", "amr", "eas", "nfe", "sas", "fin", "mid"]

FREQUENCY_COLUMNS = [
    "rsID", "variant_id", "reference_genome", "chrom", "pos",
    "ref", "alt", "overall_ac", "overall_an", "overall_af", "dataset",
] + [
    f"{pop.upper()}_{metric}"
    for pop in POPULATIONS
    for metric in ("ac", "an", "af")
]

ERROR_COLUMNS = ["rsID", "status", "error"]

QUERY = """
query VariantFrequency($rsid: String!, $dataset: DatasetId!) {
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

session = requests.Session()
session.headers.update({"Content-Type": "application/json"})


def safe_frequency(ac, an):
    if ac is None or an is None or an == 0:
        return None
    return ac / an


def fetch_variant(rsid):
    payload = {
        "query": QUERY,
        "variables": {"rsid": rsid, "dataset": DATASET},
    }

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            response = session.post(
                API_URL, json=payload, timeout=45
            )

            if response.status_code == 429 or response.status_code >= 500:
                if attempt < MAX_RETRIES:
                    retry_after = response.headers.get("Retry-After")
                    try:
                        wait = float(retry_after)
                    except (TypeError, ValueError):
                        wait = REQUEST_DELAY * attempt * 2
                    time.sleep(wait)
                    continue

            response.raise_for_status()
            result = response.json()

            if result.get("errors"):
                error_text = json.dumps(result["errors"])
                if "Multiple variants found" in error_text:
                    return None, "ambiguous_variant", error_text
                return None, "graphql_error", error_text

            variant = result.get("data", {}).get("variant")
            if variant is None:
                return None, "not_found", "No variant returned by gnomAD"

            joint = variant.get("joint") or {}
            population_data = {
                item.get("id"): item
                for item in joint.get("populations", [])
                if item.get("id")
            }

            row = {
                "rsID": rsid,
                "variant_id": variant.get("variant_id"),
                "reference_genome": variant.get("reference_genome"),
                "chrom": variant.get("chrom"),
                "pos": variant.get("pos"),
                "ref": variant.get("ref"),
                "alt": variant.get("alt"),
                "overall_ac": joint.get("ac"),
                "overall_an": joint.get("an"),
                "overall_af": safe_frequency(
                    joint.get("ac"), joint.get("an")
                ),
                "dataset": DATASET,
            }

            for pop in POPULATIONS:
                item = population_data.get(pop, {})
                ac, an = item.get("ac"), item.get("an")
                row[f"{pop.upper()}_ac"] = ac
                row[f"{pop.upper()}_an"] = an
                row[f"{pop.upper()}_af"] = safe_frequency(ac, an)

            return row, None, None

        except (requests.RequestException, ValueError) as exc:
            if attempt == MAX_RETRIES:
                return None, "request_error", str(exc)
            time.sleep(REQUEST_DELAY * attempt * 2)

    return None, "request_error", "Unknown error"


def save_outputs(rows, errors):
    pd.DataFrame(rows, columns=FREQUENCY_COLUMNS).to_csv(
        OUTPUT_FILE, index=False
    )
    pd.DataFrame(errors, columns=ERROR_COLUMNS).to_csv(
        ERROR_FILE, index=False
    )


def main():
    rsids = list(dict.fromkeys(
        line.strip()
        for line in RSID_FILE.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ))

    rows_by_id = {}
    if OUTPUT_FILE.exists():
        existing = pd.read_csv(
            OUTPUT_FILE, dtype={"rsID": str}
        )
        if "rsID" in existing.columns:
            for row in existing.to_dict("records"):
                rows_by_id[str(row["rsID"])] = row

    errors_by_id = {}
    if ERROR_FILE.exists():
        old_errors = pd.read_csv(
            ERROR_FILE, dtype=str, keep_default_na=False
        )
        if {"rsID", "status", "error"}.issubset(old_errors.columns):
            for row in old_errors.to_dict("records"):
                errors_by_id[row["rsID"]] = row

    # Existing successful results are not requested again.
    completed = set(rows_by_id)
    todo = [rsid for rsid in rsids if rsid not in completed]

    print(f"Total input rsIDs: {len(rsids)}")
    print(f"Already successful: {len(completed)}")
    print(f"To process: {len(todo)}")

    for index, rsid in enumerate(todo, start=1):
        row, status, error = fetch_variant(rsid)

        if row is not None:
            rows_by_id[rsid] = row
            errors_by_id.pop(rsid, None)
        else:
            errors_by_id[rsid] = {
                "rsID": rsid,
                "status": status,
                "error": error,
            }
            print(f"[{index}/{len(todo)}] {rsid}: {status}")

        if index % SAVE_EVERY == 0:
            save_outputs(
                list(rows_by_id.values()),
                list(errors_by_id.values()),
            )
            print(f"Saved progress: {index}/{len(todo)}")

        time.sleep(REQUEST_DELAY)

    save_outputs(
        list(rows_by_id.values()),
        list(errors_by_id.values()),
    )

    print("\nFinished.")
    print(f"Successful rsIDs: {len(rows_by_id)}")
    print(f"Unresolved/error rsIDs: {len(errors_by_id)}")
    print("\nError categories:")
    if errors_by_id:
        print(pd.Series([
            item["status"] for item in errors_by_id.values()
        ]).value_counts().to_string())
    else:
        print("None")

    print(f"\nFrequency file: {OUTPUT_FILE}")
    print(f"Error file: {ERROR_FILE}")


if __name__ == "__main__":
    main()
