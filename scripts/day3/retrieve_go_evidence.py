from Bio import Entrez
import pandas as pd
import xml.etree.ElementTree as ET
import time
from pathlib import Path

Entrez.email = input("Enter your email: ").strip()
Entrez.tool = "LS604_Day3_Practical"

data_dir = Path("data/day3")

raw_dir = data_dir / "ncbi_full_xml"
raw_dir.mkdir(parents=True, exist_ok=True)

genes = pd.read_csv(
    data_dir / "ncbi_gene_information.csv",
    dtype={"NCBI Gene ID": str}
)

records = []
status = []

for gene_id in genes["NCBI Gene ID"].unique():

    print(f"Retrieving GO evidence for {gene_id}...")

    with Entrez.efetch(
        db="gene",
        id=gene_id,
        retmode="xml"
    ) as handle:
        raw = handle.read()

    if isinstance(raw, str):
        raw = raw.encode("utf-8")

    raw_file = raw_dir / f"gene_{gene_id}.xml"
    raw_file.write_bytes(raw)

    root = ET.fromstring(raw)

    gene_count = 0

    for block in root.findall(".//Gene-commentary"):

        heading = block.findtext(
            "Gene-commentary_heading", ""
        )

        if heading.lower() != "geneontology":
            continue

        for category in block.findall(
            "./Gene-commentary_comment/Gene-commentary"
        ):

            aspect = category.findtext(
                "Gene-commentary_label", ""
            )

            for term in category.findall(
                "./Gene-commentary_comment/Gene-commentary"
            ):

                name = term.findtext(
                    "Gene-commentary_text", ""
                )

                go_ids = []

                for dbtag in term.findall(".//Dbtag"):

                    db = dbtag.findtext("Dbtag_db", "")

                    if db.lower() == "go":

                        go_id = dbtag.findtext(
                            ".//Object-id_id", ""
                        )

                        if go_id:
                            go_ids.append(
                                f"GO:{int(go_id):07d}"
                            )

                evidence = []

                for prop in term.findall(".//Other-source"):

                    text = prop.findtext(
                        "Other-source_post-text", ""
                    )

                    if text:
                        evidence.append(text)

                records.append({
                    "NCBI Gene ID": gene_id,
                    "GO category": aspect,
                    "GO ID": ";".join(go_ids),
                    "GO term": name,
                    "Evidence details": (
                        ";".join(evidence)
                        if evidence
                        else "Not extracted"
                    )
                })

                gene_count += 1

    status.append({
        "NCBI Gene ID": gene_id,
        "GO records extracted": gene_count
    })

    print("GO records extracted:", gene_count)

    time.sleep(0.5)

result = pd.DataFrame(
    records,
    columns=[
        "NCBI Gene ID",
        "GO category",
        "GO ID",
        "GO term",
        "Evidence details"
    ]
)

result.to_csv(
    data_dir / "ncbi_go_evidence.csv",
    index=False
)

status_df = pd.DataFrame(status)

status_df.to_csv(
    data_dir / "go_retrieval_status.csv",
    index=False
)

print("\nGO EVIDENCE RETRIEVAL COMPLETED")
print("-" * 50)

print(status_df.to_string(index=False))

print("\nTotal GO records:", len(result))

print("\nGO EVIDENCE PREVIEW")

if not result.empty:
    print(result.head(30).to_string(index=False))
else:
    print("No GO records extracted. Inspect raw XML.")

print("\nFiles saved successfully.")