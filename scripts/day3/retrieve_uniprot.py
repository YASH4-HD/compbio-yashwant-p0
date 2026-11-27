import requests
import pandas as pd
from pathlib import Path
from io import StringIO
from datetime import datetime, timezone

# Output directory
output_dir = Path("data/day3")
output_dir.mkdir(parents=True, exist_ok=True)

# UniProt REST API
url = "https://rest.uniprot.org/uniprotkb/stream"

# Honeybee serine protease keyword search
query = "(organism_id:7460) AND (go:0004252 OR go:0008236 OR keyword:KW-0720)"

params = {
    "query": query,
    "format": "tsv",
    "fields": (
        "accession,reviewed,protein_name,"
        "gene_names,organism_name,length,"
        "xref_interpro,xref_geneid"
    )
}

print("Retrieving UniProt protein data...")

response = requests.get(url, params=params, timeout=120)
response.raise_for_status()

df = pd.read_csv(StringIO(response.text), sep="\t")

# Save dataset
output_file = output_dir / "route1_annotation_corrected.csv"
df.to_csv(output_file, index=False)

# Display results
print("\nRETRIEVAL SUCCESSFUL")
print("-" * 45)
print("Organism: Apis mellifera")
print("Total proteins:", len(df))

print("\nReview status:")
print(df["Reviewed"].value_counts())

reviewed = df[df["Reviewed"] == "reviewed"]

print("\nReviewed proteins:")
print(reviewed[["Entry", "Protein names"]].to_string(index=False))

print("\nDataset saved to:", output_file)

# Save retrieval provenance
provenance = (
    f"Database: UniProt\n"
    f"Query: {query}\n"
    f"URL: {response.url}\n"
    f"Retrieved UTC: {datetime.now(timezone.utc).isoformat()}\n"
    f"UniProt release: {response.headers.get('X-UniProt-Release', 'Not reported')}\n"
    f"Total proteins: {len(df)}\n"
)

(output_dir / "route1_provenance.txt").write_text(
    provenance, encoding="utf-8"
)

print("Provenance record saved.")