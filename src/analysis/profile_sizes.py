#!/usr/bin/env python3
"""Size of the disease profile each individual is ranked against, per knowledge base.

For every (disease, source publication) the number of phenotype annotation lines (aspect P)
of the disease that remain after leave-one-publication-out, next to the sizes of the full
released and augmented profiles. Reads the hpoadj adjustment summaries and the corpus files.

    python profile_sizes.py            # writes _work/tables/tableS7_profile_sizes.csv
"""

import csv
from collections import defaultdict
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parent
WORK = ANALYSIS_DIR / "_work"
OUTPUT = WORK / "tables/tableS7_profile_sizes.csv"

DISEASES = {
    "OMIM:619375": "SOCS1",
    "OMIM:615513": "APDS1",
    "OMIM:616576": "NFKB1",
}

CORPORA = {
    "released_v2024": WORK / "data" / "old" / "phenotype.hpoa",
    "released_v2026": WORK / "data" / "new" / "phenotype.hpoa",
    "augmented_v2026": WORK / "hpoa_adjusted" / "new" / "phenotype_augmented.hpoa",
}

LOO_SUMMARIES = {
    "released_v2024_loo": WORK / "hpoa_adjusted" / "old" / "adjustment_summary.tsv",
    "released_v2026_loo": WORK / "hpoa_adjusted" / "new_stock" / "adjustment_summary.tsv",
    "augmented_v2026_loo": WORK / "hpoa_adjusted" / "new" / "adjustment_summary.tsv",
}


def phenotype_lines(hpoa: Path) -> dict:
    counts = defaultdict(int)
    with open(hpoa) as handle:
        for line in handle:
            if line.startswith("#") or line.startswith("database_id"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] in DISEASES and fields[10] == "P":
                counts[fields[0]] += 1
    return counts


def loo_profiles(summary: Path) -> dict:
    """(disease, pmid) -> (individuals from that publication, phenotype lines left)."""
    profiles = {}
    individuals = defaultdict(int)
    with open(summary) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["disease_id"], row["pmid"])
            individuals[key] += 1
            profiles[key] = int(row["remaining_phenotype_lines"])
    return {key: (individuals[key], profiles[key]) for key in profiles}


def main():
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    full = {name: phenotype_lines(path) for name, path in CORPORA.items()}
    loo = {name: loo_profiles(path) for name, path in LOO_SUMMARIES.items()}
    keys = sorted({key for table in loo.values() for key in table}, key=lambda k: (DISEASES[k[0]], k[1]))

    columns = ["cohort", "disease_id", "pmid_left_out", "n_individuals"] + list(CORPORA) + list(LOO_SUMMARIES)
    with open(OUTPUT, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for disease, pmid in keys:
            row = {
                "cohort": DISEASES[disease],
                "disease_id": disease,
                "pmid_left_out": pmid,
                "n_individuals": loo["augmented_v2026_loo"][(disease, pmid)][0],
            }
            for name, counts in full.items():
                row[name] = counts.get(disease, 0)
            for name, table in loo.items():
                entry = table.get((disease, pmid))
                row[name] = entry[1] if entry else full[name.removesuffix("_loo")].get(disease, 0)
            writer.writerow(row)
    print(f"wrote {OUTPUT} ({len(keys)} rows)")


if __name__ == "__main__":
    main()
