from Bio import Entrez
import pandas as pd
import time
from pathlib import Path
from datetime import datetime, timezone

Entrez.email = input("Enter your email: ").strip()
Entrez.tool = "LS604_Day3_Practical"

data_dir = Path("data/day3")

raw_dir = data_dir / "ncbi_raw"
raw_dir.mkdir(parents=True, exist_ok=True)

xml_dir = data_dir / "ncbi_full_xml"
xml_dir.mkdir(parents=True, exist_ok=True)

existing = pd.read_csv(
    data_dir / "ncbi_gene_information.csv",
    dtype={"NCBI Gene ID": str}
)

additional = pd.read_csv(
    data_dir / "additional_candidates_verified.csv",
    dtype={"NCBI Gene ID": str}
)

gene_ids = sorted(
    set(additional["NCBI Gene ID"].dropna())
    - set(existing["NCBI Gene ID"])
)

records = []

print("New genes to retrieve:", len(gene_ids))

for gene_id in gene_ids:

    print(f"\nRetrieving NCBI Gene {gene_id}")

    try:

        with Entrez.esummary(
            db="gene",
            id=gene_id,
            retmode="xml"
        ) as handle:

            raw = handle.read()

        if isinstance(raw, str):
            raw = raw.encode("utf-8")

        (raw_dir / f"gene_{gene_id}.xml").write_bytes(raw)

        with open(
            raw_dir / f"gene_{gene_id}.xml", "rb"
        ) as handle:

            data = Entrez.read(handle)

        doc = data["DocumentSummarySet"][
            "DocumentSummary"
        ][0]

        record = {
            "NCBI Gene ID": gene_id,
            "Gene symbol": str(doc.get("Name", "")),
            "Gene description": str(
                doc.get("Description", "")
            ),
            "Gene summary": str(
                doc.get("Summary", "")
                or "No summary provided"
            ),
            "Organism": str(
                doc.get("Organism", {}).get(
                    "ScientificName", ""
                )
            ),
            "NCBI Gene URL": (
                f"https://www.ncbi.nlm.nih.gov/gene/{gene_id}"
            ),
            "Retrieved UTC": datetime.now(
                timezone.utc
            ).isoformat()
        }

        time.sleep(0.5)

        with Entrez.efetch(
            db="gene",
            id=gene_id,
            retmode="xml"
        ) as handle:

            full_xml = handle.read()

        if isinstance(full_xml, str):
            full_xml = full_xml.encode("utf-8")

        (xml_dir / f"gene_{gene_id}.xml").write_bytes(
            full_xml
        )

        records.append(record)

        print("  Symbol:", record["Gene symbol"])
        print("  Description:", record["Gene description"])
        print("  Summary:", record["Gene summary"])

    except Exception as error:

        print("  Error:", error)

    time.sleep(1)

new_genes = pd.DataFrame(
    records,
    columns=existing.columns
)

new_genes.to_csv(
    data_dir / "additional_gene_information.csv",
    index=False
)

combined = pd.concat(
    [existing, new_genes],
    ignore_index=True
).drop_duplicates(
    subset=["NCBI Gene ID"]
)

combined.to_csv(
    data_dir / "ncbi_gene_information_expanded.csv",
    index=False
)

print("\nNCBI GENE RETRIEVAL COMPLETED")
print("-" * 60)

print("Additional genes retrieved:", len(new_genes))
print("Total genes in expanded dataset:", len(combined))

print(
    combined[
        [
            "NCBI Gene ID",
            "Gene symbol",
            "Gene description"
        ]
    ].to_string(index=False)
)

print("\nFiles saved successfully.")