import requests
import pandas as pd
import time
import json
from pathlib import Path
from datetime import datetime, timezone

data_dir = Path("data/day3")

raw_dir = data_dir / "interpro_raw"
raw_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(
    data_dir / "verified_candidates.csv"
)

records = []

for _, row in df.iterrows():

    accession = row["Accession"]
    length = int(row["Length"])

    print(f"Checking InterPro domain: {accession}")

    url = (
        "https://www.ebi.ac.uk/interpro/api/"
        "entry/interpro/IPR001254/"
        f"protein/uniprot/{accession}/"
    )

    try:

        response = requests.get(
            url,
            headers={"Accept": "application/json"},
            timeout=90
        )

        response.raise_for_status()

        data = response.json()

        (raw_dir / f"{accession}.json").write_text(
            json.dumps(data, indent=2),
            encoding="utf-8"
        )

        domain_count = 0

        for result in data.get("results", []):

            for protein in result.get("proteins", []):

                for location in protein.get(
                    "entry_protein_locations", []
                ):

                    fragments = location.get(
                        "fragments", []
                    )

                    if not fragments:
                        continue

                    domain_count += 1

                    positions = set()

                    for fragment in fragments:

                        start = int(fragment["start"])
                        end = int(fragment["end"])

                        positions.update(
                            range(start, end + 1)
                        )

                    domain_start = min(positions)
                    domain_end = max(positions)
                    domain_length = len(positions)

                    coverage = round(
                        domain_length / length * 100,
                        2
                    )

                    records.append({
                        "Accession": accession,
                        "Protein length": length,
                        "InterPro ID": "IPR001254",
                        "Domain number": domain_count,
                        "Domain start": domain_start,
                        "Domain end": domain_end,
                        "Domain length": domain_length,
                        "Protein coverage (%)": coverage,
                        "Retrieval status": "Retrieved",
                        "Source URL": url,
                        "Retrieved UTC": datetime.now(
                            timezone.utc
                        ).isoformat()
                    })

                    print(
                        f"  Domain {domain_count}: "
                        f"{domain_start}-{domain_end}, "
                        f"length={domain_length} aa"
                    )

        if domain_count == 0:

            records.append({
                "Accession": accession,
                "Protein length": length,
                "InterPro ID": "IPR001254",
                "Retrieval status": "No coordinates extracted",
                "Source URL": url
            })

            print("  No coordinates extracted")

    except Exception as error:

        print(f"  Retrieval error: {error}")

        records.append({
            "Accession": accession,
            "Protein length": length,
            "InterPro ID": "IPR001254",
            "Retrieval status": f"Error: {error}",
            "Source URL": url
        })

    time.sleep(1)

result = pd.DataFrame(records)

result.to_csv(
    data_dir / "interpro_domain_coordinates.csv",
    index=False
)

print("\nINTERPRO VERIFICATION COMPLETED")
print("-" * 60)

print(result.to_string(index=False))

print("\nRetrieval status:")
print(result["Retrieval status"].value_counts())

print("\nDomain coordinate file saved.")