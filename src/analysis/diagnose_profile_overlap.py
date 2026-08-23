#!/usr/bin/env python3
"""Diagnose the arm-matched (--augment-arms both) LIRICAL comparison.

Hypothesis under test: leave-one-publication-out removes thinly-supported
SPECIFIC terms from the disease profile far more often than it removes pooled
GENERAL terms. A curated patient then carries workshop terms that the LOO
profile no longer contains (only ancestor-level partial matches remain), while
the aged patient and the aged profile collapse onto the same ancestors and
match verbatim. If true, the vocabulary-only comparison is biased AGAINST the
new release, and the near-null / negative result is a floor, not the effect.

For every staged phenopacket and both arms, this reports the share of the
patient's observed terms found verbatim in (a) the full augmented profile of
the causal disease and (b) the LOO profile the patient is actually scored
against, plus the share lost to LOO.

Requires the outputs of `run_lirical_adjusted.py --augment-arms both`
(the old arm has no phenotype_augmented.hpoa after a default-mode run).

Usage, from src/analysis:
    python diagnose_profile_overlap.py
    python diagnose_profile_overlap.py --hpoa-dir _work/hpoa_adjusted.both
    python diagnose_profile_overlap.py --per-patient
"""

import argparse
import csv
import json
import re
from pathlib import Path

ANALYSIS = Path(__file__).resolve().parent
WORK = ANALYSIS / "_work"
DISEASE = {"SOCS1": "OMIM:619375", "APDS1": "OMIM:615513", "NFKB1": "OMIM:616576"}
ARMS = ("old", "new")


def hpoa_terms(path: Path, disease: str) -> set[str]:
    """Observed (non-NOT) P-aspect terms annotated to `disease` in one HPOA."""
    terms: set[str] = set()
    with open(path, encoding="utf-8") as fh:
        rows = (ln for ln in fh if not ln.startswith("#"))
        for r in csv.DictReader(rows, delimiter="\t"):
            if (r.get("database_id") == disease and r.get("aspect") == "P"
                    and (r.get("qualifier") or "").strip().upper() != "NOT"):
                terms.add(r["hpo_id"])
    return terms


def pmid_of(path: Path) -> str | None:
    try:
        meta = json.loads(path.read_text()).get("metaData", {})
        for ref in meta.get("externalReferences", []):
            if ref.get("id", "").startswith("PMID:"):
                return ref["id"].split(":", 1)[1]
    except (json.JSONDecodeError, OSError):
        pass
    m = re.match(r"PMID_(\d+)_", path.name)
    return m.group(1) if m else None


def observed_terms(path: Path) -> set[str]:
    pp = json.loads(path.read_text())
    return {f["type"]["id"] for f in pp.get("phenotypicFeatures", [])
            if not f.get("excluded", False)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hpoa-dir", type=Path, default=WORK / "hpoa_adjusted",
                    help="hpoa_adjusted directory of a --augment-arms both run")
    ap.add_argument("--per-patient", action="store_true")
    args = ap.parse_args()

    hdr = (f"{'cohort':6s} {'arm':4s} {'n':>3s} {'obs/pat':>7s} "
           f"{'in_full':>8s} {'in_loo':>7s} {'lost_by_loo':>11s}")
    print(hdr)
    print("-" * len(hdr))
    for cohort, disease in DISEASE.items():
        for arm in ARMS:
            aug = args.hpoa_dir / arm / "phenotype_augmented.hpoa"
            if not aug.exists():
                print(f"{cohort:6s} {arm:4s}   [skip] {aug} missing - "
                      "run --augment-arms both first (or point --hpoa-dir at it)")
                continue
            full = hpoa_terms(aug, disease)
            rows = []
            for pk in sorted((WORK / "lirical" / cohort / arm).glob("*.json")):
                obs = observed_terms(pk)
                if not obs:
                    continue
                pmid = pmid_of(pk)
                loo_f = args.hpoa_dir / arm / f"phenotype_PMID_{pmid}.hpoa"
                loo = hpoa_terms(loo_f, disease) if pmid and loo_f.exists() else full
                in_full = len(obs & full) / len(obs)
                in_loo = len(obs & loo) / len(obs)
                rows.append((pk.name, len(obs), in_full, in_loo, in_full - in_loo))
                if args.per_patient:
                    print(f"    {cohort} {arm} {pk.name:44s} obs={len(obs):3d} "
                          f"full={in_full:6.1%} loo={in_loo:6.1%}")
            if not rows:
                print(f"{cohort:6s} {arm:4s}   [skip] no staged phenopackets")
                continue
            n = len(rows)
            mean = lambda i: sum(r[i] for r in rows) / n  # noqa: E731
            print(f"{cohort:6s} {arm:4s} {n:3d} {mean(1):7.1f} "
                  f"{mean(2):8.1%} {mean(3):7.1%} {mean(4):11.1%}")
    print()
    print("obs/pat     = mean observed (non-excluded) terms per patient")
    print("in_full     = share of those terms present VERBATIM in the full augmented")
    print("              profile of the causal disease (before LOO)")
    print("in_loo      = same, in the LOO profile the patient is actually scored against")
    print("lost_by_loo = drop between the two, i.e. terms deleted by removing the")
    print("              patient's own publication")
    print()
    print("Reading: if in_loo is clearly lower for arm=new than arm=old, the reversal")
    print("in the arm-matched LIRICAL result is the LOO x specificity interaction")
    print("(specific terms rest on fewer publications), not a vocabulary failure.")


if __name__ == "__main__":
    main()
