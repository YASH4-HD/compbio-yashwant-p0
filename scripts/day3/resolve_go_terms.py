import requests
import pandas as pd
import time
from pathlib import Path
from datetime import datetime, timezone

data_dir = Path("data/day3")

df = pd.read_csv(
    data_dir / "ncbi_go_evidence.csv",
    dtype={"NCBI Gene ID": str}
)

go_ids = sorted(
    df["GO ID"].dropna().unique()
)

go_names = {}

for go_id in go_ids:

    print(f"Resolving {go_id}...")

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

    data = response.json()

    results = data.get("results", [])

    if results:
        name = results[0].get("name", "")
    else:
        name = "Not resolved"

    go_names[go_id] = name

    print(f"  {go_id}: {name}")

    time.sleep(0.5)

df["GO term"] = df["GO ID"].map(go_names)

df["Evidence code"] = (
    df["Evidence details"]
    .str.extract(r"evidence:\s*([A-Z]+)")
)

df.to_csv(
    data_dir / "ncbi_go_evidence_corrected.csv",
    index=False
)

mapping = pd.DataFrame([
    {
        "GO ID": go_id,
        "GO term": name,
        "Source": "QuickGO",
        "Retrieved UTC": datetime.now(
            timezone.utc
        ).isoformat()
    }
    for go_id, name in go_names.items()
])

mapping.to_csv(
    data_dir / "go_term_mapping.csv",
    index=False
)

print("\nGO TERM CORRECTION COMPLETED")
print("-" * 50)

print("Unique GO terms:", len(go_names))
print("Total annotation records:", len(df))

print("\nEvidence code distribution:")
print(df["Evidence code"].value_counts())

print("\nCorrected GO annotations:")
print(
    df[
        [
            "NCBI Gene ID",
            "GO ID",
            "GO term",
            "Evidence code"
        ]
    ].to_string(index=False)
)

print("\nCorrected files saved successfully.")