
import pandas as pd
from pathlib import Path

DATA = Path("data/day3")
RESULTS = Path("results/day3")
RESULTS.mkdir(parents=True, exist_ok=True)

# Load existing datasets
original = pd.read_csv(DATA / "verified_candidates.csv")
additional = pd.read_csv(DATA / "additional_candidates_verified.csv")
domains = pd.read_csv(DATA / "interpro_domain_coordinates_corrected.csv")
genes = pd.read_csv(DATA / "ncbi_gene_information_expanded.csv")
go = pd.read_csv(DATA / "ncbi_go_evidence_expanded.csv")
expression = pd.read_csv(DATA / "ncbi_expression_expanded.csv")
route2 = pd.read_csv(DATA / "route2_interpro.csv")

# Final selected proteins
selected = [
    "A0A7M7IJG6",
    "A0A7M7LM75",
    "A0A7M7MSZ3",
    "A0A7M7SQ79",
    "A0A7M7FY83",
    "A0A7M7GQT4",
    "A0A7M7L6X9",
    "A0A7M7FXW4",
    "A0A7M7FZ02",
    "A0A7M7GQ22",
    "A0A7M7MPA5",
    "A0A7M7ITA0"
]

# Combine protein verification records
proteins = pd.concat(
    [original, additional],
    ignore_index=True
)

proteins = proteins[
    proteins["Accession"].isin(selected)
].copy()

assert len(proteins) == 12
assert proteins["Accession"].is_unique

# Standardize gene identifiers
for df in [proteins, genes, go, expression]:
    df["NCBI Gene ID"] = (
        df["NCBI Gene ID"]
        .astype(str)
        .str.replace(r"\.0$", "", regex=True)
    )

records = []

for accession in selected:

    p = proteins[
        proteins["Accession"] == accession
    ].iloc[0]

    gene_id = p["NCBI Gene ID"]

    u = route2[
        route2["Entry"] == accession
    ].iloc[0]

    g = genes[
        genes["NCBI Gene ID"] == gene_id
    ].iloc[0]

    e = expression[
        expression["NCBI Gene ID"] == gene_id
    ].iloc[0]

    go_rows = go[
        go["NCBI Gene ID"] == gene_id
    ]

    go_terms = "; ".join(
        f"{r['GO ID']}: {r['GO term']} "
        f"({r['Evidence code']})"
        for _, r in go_rows.iterrows()
    )

    d = domains[
        domains["Accession"] == accession
    ]

    if not d.empty:
        coordinates = "; ".join(
            f"{int(r['Domain start'])}-"
            f"{int(r['Domain end'])}"
            for _, r in d.iterrows()
        )
        domain_count = d["Domain number"].nunique()

    else:
        coordinates = p["Domain coordinates"]
        domain_count = int(p["Domain count"])

    records.append({
        "UniProt accession": accession,
        "UniProt name": u["Protein names"],
        "Protein length": int(p["Length"]),
        "NCBI Gene ID": gene_id,
        "Gene symbol": g["Gene symbol"],
        "NCBI Gene description": g["Gene description"],
        "NCBI Gene summary": g["Gene summary"],
        "InterPro ID": "IPR001254",
        "Domain coordinates": coordinates,
        "Domain count": domain_count,
        "Annotated active sites": int(
            p["Annotated active sites"]
        ),
        "NCBI GO terms": go_terms,
        "GO record count": len(go_rows),
        "Expression status": e["Expression status"]
    })

master = pd.DataFrame(records)

assert len(master) == 12
assert master["UniProt accession"].is_unique
assert master["NCBI Gene ID"].nunique() == 11
assert master["Domain coordinates"].notna().all()

output = RESULTS / "master_protein_evidence.csv"
master.to_csv(output, index=False)

print("\nEVIDENCE CONSOLIDATION COMPLETED")
print("-" * 65)

print("Selected proteins:", len(master))
print("Unique genes:", master["NCBI Gene ID"].nunique())
print("Total domain matches:", master["Domain count"].sum())
print("Total annotated active sites:",
      master["Annotated active sites"].sum())

print("\nFINAL SELECTED PROTEINS")

print(master[[
    "UniProt accession",
    "NCBI Gene ID",
    "NCBI Gene description",
    "Domain coordinates",
    "Annotated active sites"
]].to_string(index=False))

print("\nSaved:", output)
