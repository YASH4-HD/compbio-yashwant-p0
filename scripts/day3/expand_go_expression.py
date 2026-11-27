import pandas as pd
import requests
import time
import re
import xml.etree.ElementTree as ET
from pathlib import Path

data_dir = Path("data/day3")
raw_dir = data_dir / "ncbi_full_xml"

genes = pd.read_csv(
    data_dir / "ncbi_gene_information_expanded.csv",
    dtype={"NCBI Gene ID": str}
)

go_records = []
expression_records = []
status_records = []

for _, row in genes.iterrows():

    gene_id = row["NCBI Gene ID"]
    symbol = row["Gene symbol"]

    print(f"\nChecking {gene_id} ({symbol})")

    raw_file = raw_dir / f"gene_{gene_id}.xml"

    if not raw_file.exists():
        print("  XML file missing")
        status_records.append({
            "NCBI Gene ID": gene_id,
            "Status": "XML missing"
        })
        continue

    root = ET.parse(raw_file).getroot()

    gene_go = []
    expression_hits = []

    for block in root.findall(".//Gene-commentary"):

        heading = block.findtext(
            "Gene-commentary_heading", ""
        )

        label = block.findtext(
            "Gene-commentary_label", ""
        )

        if heading.strip().lower() == "geneontology":

            for category in block.findall(
                "./Gene-commentary_comment/Gene-commentary"
            ):

                category_name = category.findtext(
                    "Gene-commentary_label", ""
                )

                for term in category.findall(
                    "./Gene-commentary_comment/Gene-commentary"
                ):

                    term_name = term.findtext(
                        "Gene-commentary_text", ""
                    )

                    evidence_details = "; ".join(
                        x.text.strip()
                        for x in term.findall(
                            ".//Other-source_post-text"
                        )
                        if x.text
                    )

                    evidence_codes = re.findall(
                        r"evidence:\s*([A-Z]+)",
                        evidence_details,
                        flags=re.IGNORECASE
                    )

                    for dbtag in term.findall(".//Dbtag"):

                        database = dbtag.findtext(
                            "Dbtag_db", ""
                        )

                        if database.upper() != "GO":
                            continue

                        go_number = dbtag.findtext(
                            ".//Object-id_id", ""
                        )

                        if not go_number:
                            go_number = dbtag.findtext(
                                ".//Object-id_str", ""
                            )

                        if not go_number:
                            continue

                        go_id = (
                            "GO:"
                            + go_number.replace(
                                "GO:", ""
                            ).zfill(7)
                        )

                        gene_go.append({
                            "NCBI Gene ID": gene_id,
                            "Gene symbol": symbol,
                            "GO category": category_name,
                            "GO ID": go_id,
                            "GO term": term_name,
                            "Evidence details": evidence_details,
                            "Evidence code": ";".join(
                                sorted(set(evidence_codes))
                            )
                        })

        combined = f"{heading} {label}".lower()

        keywords = [
            "expression",
            "tissue specificity",
            "tissue-specific",
            "rna-seq",
            "transcriptome"
        ]

        if any(k in combined for k in keywords):

            content = " ".join(
                x.strip()
                for x in block.itertext()
                if x.strip()
            )

            if content:
                expression_hits.append(content)

    go_records.extend(gene_go)

    expression_hits = list(
        dict.fromkeys(expression_hits)
    )

    expression_status = (
        "Expression-related information found"
        if expression_hits
        else "Not found in retrieved Gene XML"
    )

    expression_records.append({
        "NCBI Gene ID": gene_id,
        "Gene symbol": symbol,
        "Expression status": expression_status,
        "Expression information": " | ".join(
            expression_hits
        )
    })

    status_records.append({
        "NCBI Gene ID": gene_id,
        "Status": "XML examined",
        "GO records": len(gene_go)
    })

    print("  GO records:", len(gene_go))
    print("  Expression:", expression_status)

go_df = pd.DataFrame(
    go_records,
    columns=[
        "NCBI Gene ID",
        "Gene symbol",
        "GO category",
        "GO ID",
        "GO term",
        "Evidence details",
        "Evidence code"
    ]
).drop_duplicates()

# Reuse previously resolved GO names

mapping = pd.read_csv(
    data_dir / "go_term_mapping.csv"
)

go_names = dict(
    zip(mapping["GO ID"], mapping["GO term"])
)

for go_id in sorted(go_df["GO ID"].unique()):

    if go_id in go_names:
        continue

    print(f"Resolving new GO term: {go_id}")

    try:

        url = (
            "https://www.ebi.ac.uk/QuickGO/"
            f"services/ontology/go/terms/{go_id}"
        )

        response = requests.get(
            url,
            headers={"Accept": "application/json"},
            timeout=60
        )

        response.raise_for_status()

        results = response.json().get("results", [])

        if results:
            go_names[go_id] = results[0]["name"]

    except Exception as error:

        print("  Resolution error:", error)

    time.sleep(0.5)

go_df["GO term"] = (
    go_df["GO ID"].map(go_names)
    .fillna(go_df["GO term"])
)

go_df.to_csv(
    data_dir / "ncbi_go_evidence_expanded.csv",
    index=False
)

expression_df = pd.DataFrame(expression_records)

expression_df.to_csv(
    data_dir / "ncbi_expression_expanded.csv",
    index=False
)

pd.DataFrame(status_records).to_csv(
    data_dir / "expanded_retrieval_status.csv",
    index=False
)

print("\nEXPANDED EVIDENCE RETRIEVAL COMPLETED")
print("-" * 60)

print("Total genes:", len(genes))
print("Genes with GO records:",
      go_df["NCBI Gene ID"].nunique())

print("Total GO records:", len(go_df))
print("Unique GO terms:", go_df["GO ID"].nunique())

print("\nEvidence code distribution:")
print(go_df["Evidence code"].value_counts())

print("\nGO annotations:")
print(go_df[
    [
        "NCBI Gene ID",
        "GO ID",
        "GO term",
        "Evidence code"
    ]
].to_string(index=False))

print("\nExpression status:")
print(
    expression_df["Expression status"].value_counts()
)

print("\nFiles saved successfully.")