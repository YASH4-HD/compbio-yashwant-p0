import pandas as pd
import xml.etree.ElementTree as ET
from pathlib import Path

data_dir = Path("data/day3")
raw_dir = data_dir / "ncbi_full_xml"

genes = pd.read_csv(
    data_dir / "ncbi_gene_information.csv",
    dtype={"NCBI Gene ID": str}
)

records = []

for _, row in genes.iterrows():

    gene_id = row["NCBI Gene ID"]
    symbol = row["Gene symbol"]

    print(f"Checking expression: {gene_id} ({symbol})")

    raw_file = raw_dir / f"gene_{gene_id}.xml"

    if not raw_file.exists():

        records.append({
            "NCBI Gene ID": gene_id,
            "Gene symbol": symbol,
            "Expression status": "Raw XML missing",
            "Expression information": ""
        })

        continue

    root = ET.parse(raw_file).getroot()

    expression_records = []

    for block in root.findall(".//Gene-commentary"):

        heading = block.findtext(
            "Gene-commentary_heading", ""
        )

        label = block.findtext(
            "Gene-commentary_label", ""
        )

        combined = f"{heading} {label}".lower()

        keywords = [
            "expression",
            "tissue specificity",
            "tissue-specific",
            "rna-seq",
            "transcriptome"
        ]

        if any(word in combined for word in keywords):

            content = " ".join(
                text.strip()
                for text in block.itertext()
                if text.strip()
            )

            if content:
                expression_records.append(content)

    expression_records = list(
        dict.fromkeys(expression_records)
    )

    if expression_records:

        status = "Expression-related information found"
        information = " | ".join(expression_records)

    else:

        status = "Not found in retrieved Gene XML"
        information = ""

    gene_type = root.find(".//Entrezgene_type")

    gene_type_name = (
        gene_type.get("value", "")
        if gene_type is not None
        else "Not extracted"
    )

    records.append({
        "NCBI Gene ID": gene_id,
        "Gene symbol": symbol,
        "Gene type": gene_type_name,
        "Expression status": status,
        "Expression information": information,
        "NCBI Gene URL": (
            f"https://www.ncbi.nlm.nih.gov/gene/{gene_id}"
        ),
        "Raw source": str(raw_file)
    })

    print(f"  {status}")

result = pd.DataFrame(records)

result.to_csv(
    data_dir / "ncbi_expression_information.csv",
    index=False
)

print("\nEXPRESSION INFORMATION CHECK COMPLETED")
print("-" * 60)

print(
    result[
        [
            "NCBI Gene ID",
            "Gene symbol",
            "Gene type",
            "Expression status"
        ]
    ].to_string(index=False)
)

print("\nExpression status distribution:")
print(result["Expression status"].value_counts())

print("\nExpression information file saved.")