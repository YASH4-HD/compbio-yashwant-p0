
import pandas as pd
from pathlib import Path

RESULTS = Path("results/day3")
master = pd.read_csv(
    RESULTS / "master_protein_evidence.csv"
).fillna("Not provided")

# Manually collected PubMed evidence.
# No automated literature searches are performed.

evidence = {
    "A0A7M7IJG6": (
        "0 accession; 0 LOC; 1 honeybee name hit "
        "(PMID: 26529016); 5 broader name hits.",
        "Probable",
        "Annotation gap",
        "Honeybee venom proteomics reports Stubble, "
        "while UniProt retains an uncharacterized name; "
        "exact accession mapping remains provisional."
    ),

    "A0A7M7LM75": (
        "0 relevant hits in recorded searches.",
        "Probable",
        "Ambiguous",
        "UniProt names this protein proteoglycan 4, "
        "whereas NCBI identifies its gene as Filzig."
    ),

    "A0A7M7MSZ3": (
        "0 relevant hits in recorded searches.",
        "Probable",
        "Ambiguous",
        "UniProt names this protein mucin-12, "
        "whereas NCBI identifies its gene as Filzig."
    ),

    "A0A7M7SQ79": (
        "0 accession; 0 LOC; 1 honeybee pathway hit "
        "(PMID: 8563973); 4 broader protease-related "
        "hits and 7 cyclase-related hits, including "
        "PMID: 15485853.",
        "Cannot determine",
        "Ambiguous",
        "Guanylate cyclase annotations conflict with "
        "the predicted serine protease domain."
    ),

    "A0A7M7FY83": (
        "0 direct honeybee hits; 8 Drosophila "
        "Scarface search hits.",
        "Cannot determine",
        "Ambiguous",
        "The inactive protein designation and absence "
        "of annotated active sites contrast with "
        "electronically inferred protease activity."
    ),

    "A0A7M7GQT4": (
        "0 accession; 0 LOC; 3 name-search hits "
        "(PMIDs: 15050840, 17174571, 34242825).",
        "Cannot determine",
        "Ambiguous",
        "Literature supports expression or regulatory "
        "context but does not establish catalytic "
        "activity of this specific protein."
    ),

    "A0A7M7L6X9": (
        "0 hits in recorded identifier searches.",
        "Probable",
        "Genuinely unknown",
        "Two protease domains and six annotated active "
        "sites suggest possible activity, but no "
        "direct functional literature was retrieved."
    ),

    "A0A7M7FXW4": (
        "0 hits in three recorded searches.",
        "Probable",
        "Genuinely unknown",
        "A predicted protease domain and three active "
        "sites suggest activity, but no independent "
        "functional evidence was retrieved."
    ),

    "A0A7M7FZ02": (
        "0 hits in three recorded searches.",
        "Probable",
        "Genuinely unknown",
        "The modular serine protease designation is "
        "supported by computational annotations, "
        "but direct functional evidence was not found."
    ),

    "A0A7M7GQ22": (
        "0 hits in three recorded searches.",
        "Cannot determine",
        "Ambiguous",
        "The cardiomyopathy-associated protein "
        "designation and absence of annotated active "
        "sites make protease function uncertain."
    ),

    "A0A7M7MPA5": (
        "0 hits in three recorded searches.",
        "Probable",
        "Genuinely unknown",
        "A protease domain and three annotated active "
        "sites suggest activity, but no direct "
        "functional literature was retrieved."
    ),

    "A0A7M7ITA0": (
        "0 hits in three recorded searches.",
        "Probable",
        "Genuinely unknown",
        "Computational evidence suggests protease "
        "activity, but no direct functional "
        "literature was retrieved."
    )
}

# Additional evidence required for ambiguous proteins

additional_evidence = {
    "A0A7M7LM75":
        "Verify gene-to-protein mapping and perform "
        "sequence alignment and catalytic assays.",

    "A0A7M7MSZ3":
        "Compare the two Filzig protein records "
        "and validate their functional annotations.",

    "A0A7M7SQ79":
        "Examine domain architecture and test "
        "protease and cyclase activities separately.",

    "A0A7M7FY83":
        "Align catalytic residues with active "
        "proteases and perform an activity assay.",

    "A0A7M7GQT4":
        "Verify catalytic residues and perform "
        "a substrate cleavage assay.",

    "A0A7M7GQ22":
        "Validate domain architecture, catalytic "
        "residues and gene-to-protein mapping."
}

rows = []

classification_lines = [
    "# Protein classifications",
    "",
    "Organism: Apis mellifera",
    "",
    "Domain: IPR001254",
    "",
    "Classifications are provisional and limited "
    "to the retrieved evidence.",
    "",
    "Probable activity indicates computational "
    "support, not experimental validation.",
    ""
]

for _, r in master.iterrows():

    accession = r["UniProt accession"]

    hits, activity, category, reason = evidence[accession]

    summary = str(r["NCBI Gene summary"])

    gene_text = (
        f"Gene ID: {r['NCBI Gene ID']}; "
        f"Description: {r['NCBI Gene description']}; "
        f"Summary: {summary}; "
        f"Expression: {r['Expression status']}"
    )

    domain_text = (
        f"IPR001254; residues "
        f"{r['Domain coordinates']}; "
        f"{r['Domain count']} domain match(es); "
        f"{r['Annotated active sites']} "
        f"annotated active site(s)"
    )

    rows.append({
        "UniProt accession": accession,
        "UniProt name": r["UniProt name"],
        "InterPro domain hit": domain_text,
        "NCBI Gene summary": gene_text,
        "NCBI GO terms (best)": r["NCBI GO terms"],
        "PubMed hits": hits,
        "Evidence for catalytic activity": activity,
        "Classification": category
    })

    classification_lines.extend([
        f"## {accession}",
        "",
        f"**Protein:** {r['UniProt name']}",
        "",
        f"**Classification:** {category}",
        "",
        f"**Justification:** {reason}",
        ""
    ])

    if accession in additional_evidence:
        classification_lines.extend([
            "**Additional evidence required:** "
            + additional_evidence[accession],
            ""
        ])

# FILE 1: Annotation table

table = pd.DataFrame(rows)

assert table.shape == (12, 8)
assert table["UniProt accession"].is_unique

table.to_csv(
    RESULTS / "annotation_table.csv",
    index=False
)

# FILE 2: Classifications

classification_lines.extend([
    "## Limitations",
    "",
    "Zero PubMed hits do not establish that no "
    "research exists. Alternative identifiers and "
    "full-text literature may contain evidence.",
    "",
    "Two selected protein accessions map to the "
    "same Filzig gene and are not independent "
    "gene-level observations.",
    "",
    "All retrieved GO evidence codes were IEA. "
    "Annotated active sites are not proof of "
    "catalytic activity."
])

(RESULTS / "classifications.md").write_text(
    "\n".join(classification_lines),
    encoding="utf-8"
)

# FILE 3: Evidence summary

summary = """# Evidence summary

A0A7M7IJG6 — Uncharacterized protein LOC410402 isoform X1.

This Apis mellifera protein contains the IPR001254
serine protease domain between residues 493 and 735.
UniProt records three annotated active sites, suggesting
potential catalytic activity. NCBI Gene identifies the
corresponding gene, 410402, as serine proteinase stubble.
The retrieved Gene Ontology annotations include serine-type
endopeptidase activity and proteolysis, both supported by
inferred electronic annotation rather than direct
experimental evidence.

Manual PubMed searches using the accession and LOC
identifier returned zero results. However, a protein-name
search identified a honeybee venom proteomics study
(PMID: 26529016) reporting Stubble in queen venom.
Comparative studies in Drosophila also associate Stubble
with developmental processes and protease-domain function.

The difference between the uncharacterized UniProt name
and the NCBI description suggests a potential annotation
gap. Nevertheless, the exact mapping between the published
honeybee protein and this accession remains unverified,
and no direct catalytic assay was identified.

Confidence is moderate for the provisional annotation-gap
classification and predicted activity. Verified sequence
mapping, catalytic residue alignment and a biochemical
substrate cleavage assay could strengthen or change
this interpretation.
"""

(RESULTS / "evidence_summary.md").write_text(
    summary,
    encoding="utf-8"
)

# FILE 4: Methods note

methods = """# Methods note

Python, pandas, requests and Biopython were used to retrieve
and consolidate protein and gene information from UniProt,
InterPro, NCBI Gene and QuickGO.
The UniProt annotation-based search returned 162 proteins,
while the IPR001254 domain-based search returned 113 proteins,
from which 12 candidates representing 11 genes were selected.
PubMed searches were performed manually using protein
accessions, LOC identifiers and alternative protein names.
AI assistance was used to develop Python scripts and organize
evidence, while literature searches and interpretation were
reviewed manually.
An initial keyword-only UniProt query retrieved 108 proteins
and was corrected by including relevant GO terms, and an
initial InterPro parser was corrected after inspecting the
raw JSON structure.
The analysis distinguishes computational annotations from
experimental evidence and treats zero literature hits as
a search limitation rather than proof of absent function.
"""

(RESULTS / "methods_note.md").write_text(
    methods,
    encoding="utf-8"
)

# Final validation

print("\nSUBMISSION FILES GENERATED")
print("-" * 60)

for filename in [
    "annotation_table.csv",
    "classifications.md",
    "evidence_summary.md",
    "methods_note.md"
]:
    path = RESULTS / filename
    print(filename, "|", path.stat().st_size, "bytes")

print("\nAnnotation table dimensions:", table.shape)

print("\nCLASSIFICATION DISTRIBUTION")
print(table["Classification"].value_counts())

print("\nCATALYTIC ACTIVITY DISTRIBUTION")
print(table["Evidence for catalytic activity"].value_counts())

print("\nAll four submission files created successfully.")
