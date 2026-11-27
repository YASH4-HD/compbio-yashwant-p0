from Bio import Entrez
import pandas as pd
import time
import json
from pathlib import Path
from datetime import datetime, timezone
from urllib.parse import quote_plus

# NCBI requires a contact email
Entrez.email = input("Enter your email: ").strip()
Entrez.tool = "LS604_Day3_Practical"

data_dir = Path("data/day3")

raw_dir = data_dir / "ncbi_raw"
raw_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(
    data_dir / "verified_candidates.csv",
    dtype={"NCBI Gene ID": str}
)

gene_ids = sorted({
    gene.strip()
    for value in df["NCBI Gene ID"].dropna()
    for gene in value.split(";")
    if gene.strip()
})

records = []

for gene_id in gene_ids:

    print(f"Retrieving NCBI Gene {gene_id}...")

    with Entrez.esummary(
        db="gene",
        id=gene_id,
        retmode="xml"
    ) as handle:

        raw = handle.read()

    raw_file = raw_dir / f"gene_{gene_id}.xml"

    raw_file.write_bytes(
        raw if isinstance(raw, bytes)
        else raw.encode("utf-8")
    )

    with raw_file.open("rb") as handle:
        data = Entrez.read(handle)

    documents = data[
        "DocumentSummarySet"
    ]["DocumentSummary"]

    if not documents:
        print(f"No record found for {gene_id}")
        continue

    gene = documents[0]

    symbol = str(gene.get("Name", ""))
    description = str(gene.get("Description", ""))
    summary = str(gene.get("Summary", ""))

    records.append({
        "NCBI Gene ID": gene_id,
        "Gene symbol": symbol,
        "Gene description": description,
        "Gene summary": (
            summary if summary
            else "No summary provided"
        ),
        "Organism": str(
            gene.get("Organism", {}).get(
                "ScientificName", ""
            )
        ),
        "NCBI Gene URL": (
            f"https://www.ncbi.nlm.nih.gov/gene/"
            f"{gene_id}"
        ),
        "Retrieved UTC": datetime.now(
            timezone.utc
        ).isoformat()
    })

    time.sleep(0.5)

genes = pd.DataFrame(records)

genes.to_csv(
    data_dir / "ncbi_gene_information.csv",
    index=False
)

# Combine gene information with protein records
merged = df.merge(
    genes,
    on="NCBI Gene ID",
    how="left"
)

# Prepare links for MANUAL PubMed searches
def pubmed_link(row):

    accession = row["Accession"]
    gene_id = row["NCBI Gene ID"]

    query = (
        f'("{accession}" OR "LOC{gene_id}") '
        'AND "Apis mellifera"'
    )

    return (
        "https://pubmed.ncbi.nlm.nih.gov/"
        "?term=" + quote_plus(query)
    )

merged["Manual PubMed search URL"] = (
    merged.apply(pubmed_link, axis=1)
)

merged.to_csv(
    data_dir / "protein_gene_evidence.csv",
    index=False
)

print("\nNCBI GENE RETRIEVAL COMPLETED")
print("-" * 50)

print("Unique genes retrieved:", len(genes))
print("Protein records:", len(merged))

print(
    genes[
        [
            "NCBI Gene ID",
            "Gene symbol",
            "Gene description"
        ]
    ].to_string(index=False)
)

print("\nGene summaries and search links saved.")
print("PubMed searches must be performed manually.")