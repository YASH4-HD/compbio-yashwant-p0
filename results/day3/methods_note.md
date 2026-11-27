# Methods note

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
