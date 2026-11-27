import requests
import pandas as pd
import time
from pathlib import Path

data_dir = Path("data/day3")

df = pd.read_csv(
    data_dir / "preliminary_candidates.csv"
)

records = []

for _, row in df.iterrows():

    accession = row["Entry"]

    print(f"Checking {accession}...")

    url = (
        "https://rest.uniprot.org/"
        f"uniprotkb/{accession}.json"
    )

    response = requests.get(url, timeout=60)
    response.raise_for_status()

    data = response.json()

    refs = data.get(
        "uniProtKBCrossReferences", []
    )

    gene_ids = [
        ref["id"]
        for ref in refs
        if ref["database"] == "GeneID"
    ]

    interpro_hits = [
        ref["id"]
        for ref in refs
        if ref["database"] == "InterPro"
    ]

    active_sites = [
        feature
        for feature in data.get("features", [])
        if feature.get("type") == "Active site"
    ]

    records.append({
        "Accession": accession,
        "Protein name": row["Protein names"],
        "Length": row["Length"],
        "NCBI Gene ID": ";".join(gene_ids),
        "IPR001254 present": (
            "IPR001254" in interpro_hits
        ),
        "Annotated active sites": len(active_sites),
        "UniProt URL": (
            f"https://www.uniprot.org/uniprotkb/"
            f"{accession}/entry"
        )
    })

    time.sleep(0.5)

result = pd.DataFrame(records)

result.to_csv(
    data_dir / "verified_candidates.csv",
    index=False
)

print("\nCANDIDATE VERIFICATION")
print("-" * 60)

print(
    result[
        [
            "Accession",
            "NCBI Gene ID",
            "IPR001254 present",
            "Annotated active sites"
        ]
    ].to_string(index=False)
)

unique_genes = {
    gene
    for value in result["NCBI Gene ID"]
    for gene in value.split(";")
    if gene
}

print("\nTotal proteins:", len(result))
print("Unique mapped Gene IDs:", len(unique_genes))

print(
    "Proteins without Gene ID:",
    (result["NCBI Gene ID"] == "").sum()
)

print("\nVerification file saved successfully.")