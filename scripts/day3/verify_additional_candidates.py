import requests
import pandas as pd
import time
import json
from pathlib import Path

data_dir = Path("data/day3")
raw_dir = data_dir / "additional_candidates_raw"
raw_dir.mkdir(parents=True, exist_ok=True)

accessions = [
    "A0A7M7FXW4",
    "A0A7M7FZ02",
    "A0A7M7LM75",
    "A0A7M7MSZ3",
    "A0A7M7SQ79",
    "A0A7M7GQ22"
]

records = []

for accession in accessions:

    print(f"\nChecking {accession}")

    try:

        url = (
            "https://rest.uniprot.org/"
            f"uniprotkb/{accession}.json"
        )

        response = requests.get(url, timeout=60)
        response.raise_for_status()
        data = response.json()

        (raw_dir / f"{accession}_uniprot.json").write_text(
            json.dumps(data, indent=2),
            encoding="utf-8"
        )

        length = data["sequence"]["length"]

        taxon = data.get("organism", {}).get("taxonId")

        refs = data.get("uniProtKBCrossReferences", [])

        gene_ids = [
            x["id"] for x in refs
            if x["database"] == "GeneID"
        ]

        active_sites = [
            x for x in data.get("features", [])
            if x["type"] == "Active site"
        ]

        interpro_url = (
            "https://www.ebi.ac.uk/interpro/api/"
            "entry/interpro/IPR001254/"
            f"protein/uniprot/{accession}/"
        )

        response = requests.get(
            interpro_url,
            headers={"Accept": "application/json"},
            timeout=90
        )

        response.raise_for_status()
        ip_data = response.json()

        (raw_dir / f"{accession}_interpro.json").write_text(
            json.dumps(ip_data, indent=2),
            encoding="utf-8"
        )

        domains = []

        for protein in ip_data.get("proteins", []):

            if protein.get("accession", "").upper() != accession:
                continue

            for location in protein.get(
                "entry_protein_locations", []
            ):

                fragments = location.get("fragments", [])

                if not fragments:
                    continue

                positions = set()

                for fragment in fragments:

                    positions.update(range(
                        int(fragment["start"]),
                        int(fragment["end"]) + 1
                    ))

                domains.append({
                    "start": min(positions),
                    "end": max(positions),
                    "length": len(positions)
                })

        coordinates = "; ".join(
            f"{d['start']}-{d['end']}"
            for d in domains
        )

        domain_lengths = "; ".join(
            str(d["length"]) for d in domains
        )

        records.append({
            "Accession": accession,
            "Length": length,
            "Taxonomy ID": taxon,
            "NCBI Gene ID": ";".join(gene_ids),
            "Annotated active sites": len(active_sites),
            "Domain count": len(domains),
            "Domain coordinates": coordinates,
            "Domain lengths": domain_lengths,
            "Retrieval status": "Retrieved"
        })

        print("  Gene ID:", gene_ids)
        print("  Active sites:", len(active_sites))
        print("  Domain coordinates:", coordinates)
        print("  Domain lengths:", domain_lengths)

    except Exception as error:

        print("  Error:", error)

        records.append({
            "Accession": accession,
            "Retrieval status": str(error)
        })

    time.sleep(1)

result = pd.DataFrame(records)

result.to_csv(
    data_dir / "additional_candidates_verified.csv",
    index=False
)

print("\nADDITIONAL CANDIDATE VERIFICATION")
print("-" * 60)
print(result.to_string(index=False))

print("\nResults saved successfully.")