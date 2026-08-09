#!/usr/bin/env python3
"""Run the LIRICAL benchmark on HPOA files that carry the ESID cohort annotations.

The workshop contributed both new HPO terms and the curated cohorts. The 2026-06-23 HPOA
release does not contain the cohort curation yet, so the new arm is built here:

  1. `hpotools hpoadjust --augment` adds one annotation per disease and observed term,
     with the frequency pooled over all cohort publications (ancestor-aware counting) and
     every contributing PMID in the reference field
  2. for each publication, the same command writes a copy of that HPOA in which the
     publication's own contribution is subtracted again (leave-one-publication-out)
  3. one data bundle per publication is assembled (symlinks + the adjusted phenotype.hpoa)
  4. `lirical benchmark` runs per (cohort, arm, bundle) group
  5. the group CSVs are merged into _work/lirical/<cohort>.<arm>.csv, the files the
     notebook (esid4hpo_analysis.ipynb, section 4) loads unchanged

`--augment-arms` decides which arm receives the cohort annotations:

  new   (default) new arm = 2026 terms + cohort annotations, old arm = 2024 HPOA as published.
                  Measures the workshop as a whole: new terms together with the curation.
  both            both arms carry the cohort annotations, the old arm from the aged
                  phenopackets. Isolates the effect of the new terms alone.
  none            no augmentation, only the leave-one-publication-out removal.

Usage, from src/analysis:
    python run_lirical_adjusted.py
    python run_lirical_adjusted.py --dry-run
    python run_lirical_adjusted.py --augment-arms both

Requires staged phenopackets in _work/lirical/<cohort>/{old,new}/ (notebook section 3)
and the version-matched data bundles in _work/data/{old,new}/.
"""

import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parent
WORK = ANALYSIS_DIR / "_work"
COHORT_DIR = ANALYSIS_DIR.parent / "cohorts"
ARMS = ("old", "new")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--hpotools-jar", type=Path,
                        default=ANALYSIS_DIR.parents[2] / "hpotools" / "target" / "hpotools.jar")
    parser.add_argument("--lirical-jar", type=Path,
                        default=ANALYSIS_DIR / "_tools" / "lirical-cli-2.4.1" / "lirical-cli-2.4.1.jar")
    parser.add_argument("--augment-arms", choices=["new", "both", "none"], default="new")
    parser.add_argument("--cohorts", nargs="+", default=["APDS", "NFKB1", "SOCS1"])
    parser.add_argument("--skip-adjust", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def run(cmd, dry_run):
    print("  $", " ".join(str(c) for c in cmd))
    if not dry_run:
        subprocess.run([str(c) for c in cmd], check=True)


def augments(arm, choice):
    return choice == "both" or (choice == "new" and arm == "new")


def cohort_source(arm, cohorts, dry_run):
    """Curated phenopackets for the new arm, their aged counterparts for the old arm."""
    if arm == "new":
        return COHORT_DIR
    staged = WORK / "hpoa_adjusted" / "aged_cohorts"
    if dry_run:
        return staged
    if staged.exists():
        shutil.rmtree(staged)
    staged.mkdir(parents=True)
    for cohort in cohorts:
        for packet in sorted((WORK / "lirical" / cohort / "old").glob("*.json")):
            shutil.copy(packet, staged / packet.name)
    print(f"[stage:{arm}] {len(list(staged.glob('*.json')))} aged phenopackets in {staged}")
    return staged


def adjust_hpoas(args):
    for arm in ARMS:
        outdir = WORK / "hpoa_adjusted" / arm
        if args.skip_adjust and (outdir / "adjustment_summary.tsv").exists():
            print(f"[adjust:{arm}] reusing {outdir}")
            continue
        augmented = augments(arm, args.augment_arms)
        print(f"[adjust:{arm}] augment={augmented}")
        cmd = ["java", "-jar", args.hpotools_jar, "hpoadjust",
               "-a", WORK / "data" / arm / "phenotype.hpoa",
               "--hpo", WORK / "data" / arm / "hp.json",
               "-p", cohort_source(arm, args.cohorts, args.dry_run) if augmented else COHORT_DIR,
               "-o", outdir]
        if augmented:
            cmd.append("--augment")
        run(cmd, args.dry_run)


def read_summary(arm):
    path = WORK / "hpoa_adjusted" / arm / "adjustment_summary.tsv"
    adjusted, insufficient = {}, []
    if not path.exists():
        print(f"[warn] {path} missing (dry run before first adjust?)")
        return adjusted, insufficient
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["adjusted_hpoa"] != "-":
                adjusted[row["pmid"]] = Path(row["adjusted_hpoa"])
            if row["status"] == "INSUFFICIENT_DATA":
                insufficient.append((row["phenopacket_id"], row["pmid"], row["disease_id"]))
    return adjusted, insufficient


def build_bundles(arm, adjusted, dry_run):
    stock = WORK / "data" / arm
    bundles = {}
    for pmid, hpoa in adjusted.items():
        bundle = WORK / "hpoa_adjusted" / arm / ("bundle_" + pmid.replace(":", "_"))
        bundles[pmid] = bundle
        if dry_run:
            print(f"[bundle:{arm}] would build {bundle}")
            continue
        if bundle.exists():
            shutil.rmtree(bundle)
        bundle.mkdir(parents=True)
        for entry in stock.iterdir():
            if entry.name != "phenotype.hpoa":
                (bundle / entry.name).symlink_to(entry.resolve())
        shutil.copy(hpoa, bundle / "phenotype.hpoa")
        print(f"[bundle:{arm}] {bundle.name} <- {hpoa.name}")
    return bundles


def phenopacket_pmid(path):
    try:
        meta = json.loads(path.read_text()).get("metaData", {})
        for ref in meta.get("externalReferences", []):
            if ref.get("id", "").startswith("PMID:"):
                return ref["id"]
    except (json.JSONDecodeError, OSError):
        pass
    match = re.match(r"PMID_(\d+)_", path.name)
    return f"PMID:{match.group(1)}" if match else None


def benchmark_cohort(cohort, arm, bundles, args):
    staged = WORK / "lirical" / cohort / arm
    packets = sorted(staged.glob("*.json"))
    if not packets:
        print(f"[lirical:{cohort}:{arm}] no staged phenopackets in {staged}, skipping")
        return
    groups = defaultdict(list)
    for packet in packets:
        pmid = phenopacket_pmid(packet)
        groups[pmid if pmid in bundles else "stock"].append(packet)

    parts_dir = WORK / "lirical" / "parts"
    parts_dir.mkdir(parents=True, exist_ok=True)
    parts = []
    for label, members in sorted(groups.items(), key=lambda item: item[0]):
        data_dir = WORK / "data" / arm if label == "stock" else bundles[label]
        part = parts_dir / f"{cohort}.{arm}.{label.replace(':', '_')}.csv"
        print(f"[lirical:{cohort}:{arm}] {label}: {len(members)} phenopackets vs {data_dir.name}")
        run(["java", "-jar", args.lirical_jar, "benchmark",
             "-d", data_dir, "-o", part] + members, args.dry_run)
        parts.append(part)

    merged = WORK / "lirical" / f"{cohort}.{arm}.csv"
    if args.dry_run:
        print(f"[merge] would write {merged} from {len(parts)} part(s)")
        return
    with open(merged, "w", newline="") as out:
        writer = None
        for part in parts:
            with open(part, newline="") as handle:
                reader = csv.reader(handle)
                header = next(reader)
                if writer is None:
                    writer = csv.writer(out)
                    writer.writerow(header)
                writer.writerows(reader)
    print(f"[merge] wrote {merged}")


def main():
    args = parse_args()
    for jar in (args.hpotools_jar, args.lirical_jar):
        if not args.dry_run and not jar.exists():
            sys.exit(f"jar not found: {jar}")

    adjust_hpoas(args)
    for arm in ARMS:
        adjusted, insufficient = read_summary(arm)
        bundles = build_bundles(arm, adjusted, args.dry_run)
        for cohort in args.cohorts:
            benchmark_cohort(cohort, arm, bundles, args)
        if insufficient:
            seen = sorted({(pmid, disease) for _, pmid, disease in insufficient})
            print(f"[note:{arm}] {len(insufficient)} phenopacket(s) have no phenotype annotation left after "
                  "removing their publication: "
                  + "; ".join(f"{pmid} is the sole source of {disease}" for pmid, disease in seen))

    print("\nDone. Re-run the notebook from section 4 to rebuild rank tables and figures.")


if __name__ == "__main__":
    main()
