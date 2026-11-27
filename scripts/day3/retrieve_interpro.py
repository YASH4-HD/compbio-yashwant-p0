import requests
import pandas as pd
from pathlib import Path
from io import StringIO
from datetime import datetime, timezone

output_dir = Path("data/day3")
output_dir.mkdir(parents=True, exist_ok=True)

url = "https://rest.uniprot.org/uniprotkb/stream"

query = "(organism_id:7460) AND (xref:IPR001254)"

params = {
    "query": query,
    "format": "tsv",
    "fields": (
        "accession,reviewed,protein_name,"
        "gene_names,organism_name,length,"
        "xref_interpro,xref_geneid"
    )
}

print("Retrieving IPR001254 domain-matched proteins...")

response = requests.get(
    url, params=params, timeout=120
)
response.raise_for_status()

df = pd.read_csv(
    StringIO(response.text), sep="\t"
)

output_file = output_dir / "route2_interpro.csv"
df.to_csv(output_file, index=False)

print("\nRETRIEVAL SUCCESSFUL")
print("-" * 45)
print("Organism: Apis mellifera")
print("InterPro domain: IPR001254")
print("Total proteins:", len(df))

print("\nReview status:")
print(df["Reviewed"].value_counts())

reviewed = df[df["Reviewed"] == "reviewed"]

print("\nReviewed proteins:")
print(
    reviewed[["Entry", "Protein names"]]
    .to_string(index=False)
)

# Compare both retrieval routes
route1 = pd.read_csv(
    output_dir / "route1_annotation_corrected.csv"
)

a = set(route1["Entry"])
b = set(df["Entry"])

print("\nCOMPARISON OF BOTH ROUTES")
print("-" * 45)
print("Annotation route:", len(a))
print("Domain route:", len(b))
print("Common proteins:", len(a & b))
print("Annotation exclusive:", len(a - b))
print("Domain exclusive:", len(b - a))
print("Total unique proteins:", len(a | b))

# Save provenance
provenance = (
    "Database: UniProtKB InterPro cross-references\n"
    f"Query: {query}\n"
    f"URL: {response.url}\n"
    f"Retrieved UTC: {datetime.now(timezone.utc).isoformat()}\n"
    f"UniProt release: {response.headers.get('X-UniProt-Release', 'Not reported')}\n"
    f"Total proteins: {len(df)}\n"
)

(output_dir / "route2_provenance.txt").write_text(
    provenance, encoding="utf-8"
)

print("\nDataset saved to:", output_file)
print("Provenance record saved.")