import pandas as pd
from pathlib import Path

data_dir = Path("data/day3")

# Load both datasets
route1 = pd.read_csv(
    data_dir / "route1_annotation_corrected.csv"
)

route2 = pd.read_csv(
    data_dir / "route2_interpro.csv"
)

# Select unreviewed domain-matched proteins
candidates = route2[
    route2["Reviewed"] == "unreviewed"
].copy()

# Identify potentially incomplete annotations
patterns = (
    r"uncharacterized|unknown|hypothetical|"
    r"putative|predicted|inactive|"
    r"masquerade|scarface"
)

candidates["Screening flag"] = (
    candidates["Protein names"]
    .str.contains(
        patterns,
        case=False,
        na=False,
        regex=True
    )
)

flagged = candidates[
    candidates["Screening flag"]
].copy()

# Save screening results
candidates.to_csv(
    data_dir / "all_unreviewed_candidates.csv",
    index=False
)

flagged.to_csv(
    data_dir / "preliminary_candidates.csv",
    index=False
)

# Save unique reviewed reference proteins
reviewed = pd.concat([
    route1[route1["Reviewed"] == "reviewed"],
    route2[route2["Reviewed"] == "reviewed"]
]).drop_duplicates(subset="Entry")

reviewed.to_csv(
    data_dir / "reviewed_reference_proteins.csv",
    index=False
)

print("\nCANDIDATE SCREENING RESULTS")
print("-" * 50)

print("Domain-matched proteins:", len(route2))
print("Unreviewed proteins:", len(candidates))
print("Flagged candidates:", len(flagged))
print("Unique reviewed references:", len(reviewed))

print("\nPRELIMINARY CANDIDATES")
print("-" * 50)

print(
    flagged[
        ["Entry", "Protein names", "Length"]
    ].to_string(index=False)
)

print("\nScreening files saved successfully.")
print("Final classification requires evidence review.")