import json
import pandas as pd
from pathlib import Path

data_dir = Path("data/day3")
raw_dir = data_dir / "interpro_raw"

df = pd.read_csv(
    data_dir / "verified_candidates.csv"
)

records = []

for _, row in df.iterrows():

    accession = row["Accession"]
    length = int(row["Length"])

    print(f"Extracting coordinates: {accession}")

    raw_file = raw_dir / f"{accession}.json"

    if not raw_file.exists():
        print("  Raw JSON file missing")
        continue

    data = json.loads(
        raw_file.read_text(encoding="utf-8")
    )

    domain_count = 0

    for protein in data.get("proteins", []):

        if protein.get("accession", "").upper() != accession.upper():
            continue

        for location in protein.get(
            "entry_protein_locations", []
        ):

            fragments = location.get("fragments", [])

            if not fragments:
                continue

            domain_count += 1

            positions = set()
            fragment_status = []

            for fragment in fragments:

                start = int(fragment["start"])
                end = int(fragment["end"])

                positions.update(
                    range(start, end + 1)
                )

                fragment_status.append(
                    fragment.get("dc-status", "Unknown")
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
                "Fragment status": ";".join(
                    fragment_status
                ),
                "Match score": location.get("score"),
                "Retrieval status": "Retrieved",
                "Raw source": str(raw_file)
            })

            print(
                f"  Domain {domain_count}: "
                f"{domain_start}-{domain_end}, "
                f"{domain_length} aa, "
                f"coverage={coverage}%"
            )

    if domain_count == 0:

        records.append({
            "Accession": accession,
            "Protein length": length,
            "InterPro ID": "IPR001254",
            "Retrieval status": "No coordinates extracted"
        })

result = pd.DataFrame(records)

result.to_csv(
    data_dir / "interpro_domain_coordinates_corrected.csv",
    index=False
)

successful = result[
    result["Retrieval status"] == "Retrieved"
]

print("\nINTERPRO EXTRACTION COMPLETED")
print("-" * 60)

print("Candidate proteins:", len(df))

print(
    "Proteins with coordinates:",
    successful["Accession"].nunique()
)

print("Total domain matches:", len(successful))

print("\nDOMAIN COORDINATES")

if not successful.empty:

    print(
        successful[
            [
                "Accession",
                "Domain number",
                "Domain start",
                "Domain end",
                "Domain length",
                "Protein coverage (%)",
                "Fragment status"
            ]
        ].to_string(index=False)
    )

print("\nCorrected coordinate file saved.")