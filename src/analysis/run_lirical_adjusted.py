#!/usr/bin/env python3
"""Run the LIRICAL benchmark for the ESID4HPO comparisons.

The workshop contributed both new HPO terms (the vocabulary) and the curated cohorts
(the disease annotations). The 2026-06-23 HPOA release does not contain the cohort
curation yet, so the augmented knowledge base is built here with hpoadj:

  1. `hpoadj augment` adds one annotation per disease and observed term, with the
     frequency pooled over all cohort publications (ancestor-aware counting) and every
     contributing PMID in the reference field
  2. `hpoadj loo` writes, for each publication, a copy of an HPOA in which the
     publication's own contribution is subtracted again (leave-one-publication-out)
  3. one data bundle per publication is assembled (symlinks + the adjusted phenotype.hpoa)
  4. `lirical benchmark` runs per (cohort, vocabulary, knowledge base, bundle) group
  5. the group CSVs are merged into _work/lirical/<cohort>.<job>.csv, which the
     notebook (esid4hpo_analysis.ipynb, section 4) loads

Two things vary between arms: the VOCABULARY the individual is described in
(old = aged phenopackets in v2024-08-13 terms, new = curated phenopackets in
v2026-06-23 terms) and the KNOWLEDGE BASE ranked against (always leave-one-
publication-out):

  old_stock   v2024-08-13 hp.json + released phenotype.hpoa      -> _work/hpoa_adjusted/old
  new_stock   v2026-06-23 hp.json + released phenotype.hpoa      -> _work/hpoa_adjusted/new_stock
  new_aug     v2026-06-23 hp.json + phenotype.hpoa augmented     -> _work/hpoa_adjusted/new
              with the cohort annotations

Four (vocabulary, knowledge base) jobs cover the three comparisons the manuscript reports:

  job                        merged CSV                          used by
  old   x old_stock          <cohort>.old.csv                    main (pre-workshop arm)
  new   x new_aug            <cohort>.new.csv                    main, vocab, corpus (post-workshop arm)
  old   x new_aug            <cohort>.old_terms_vs_new_hpoa.csv  vocab  (Fig 4A: knowledge base fixed, vocabulary varies)
  new   x new_stock          <cohort>.new_terms_vs_stock_hpoa.csv corpus (Fig 4B: vocabulary fixed, knowledge base varies)

  main    old x old_stock  vs  new x new_aug   the ESID4HPO output as a whole (Supplementary Figure)
  vocab   old x new_aug    vs  new x new_aug   effect of the new and re-parented terms alone
  corpus  new x new_stock  vs  new x new_aug   effect of the added disease annotations alone

Usage, from src/analysis:
    python run_lirical_adjusted.py                       # all three comparisons (4 jobs)
    python run_lirical_adjusted.py --comparisons main    # the two headline jobs only
    python run_lirical_adjusted.py --comparisons vocab corpus --skip-adjust
    python run_lirical_adjusted.py --dry-run

Requires staged phenopackets in _work/lirical/<cohort>/{old,new}/ (notebook section 3)
and the version-matched data bundles in _work/data/{old,new}/.

The earlier `--augment-arms both` mode (both arms carrying their own-vocabulary cohort
annotations) is superseded by the `vocab` comparison and no longer offered; its outputs,
if present, remain under _work/hpoa_adjusted.both/ and _work/lirical/<cohort>.<arm>.both.csv.
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

# knowledge base -> (release arm of the data bundle, augmented with the cohort, output dir)
KNOWLEDGE_BASES = {
    "old_stock": dict(release="old", augment=False, outdir="old"),
    "new_stock": dict(release="new", augment=False, outdir="new_stock"),
    "new_aug":   dict(release="new", augment=True,  outdir="new"),
}

# job -> (vocabulary of the staged phenopackets, knowledge base, merged CSV name)
JOBS = {
    "old_x_old_stock": ("old", "old_stock", "{cohort}.old.csv"),
    "new_x_new_aug":   ("new", "new_aug",   "{cohort}.new.csv"),
    "old_x_new_aug":   ("old", "new_aug",   "{cohort}.old_terms_vs_new_hpoa.csv"),
    "new_x_new_stock": ("new", "new_stock", "{cohort}.new_terms_vs_stock_hpoa.csv"),
}

COMPARISONS = {
    "main":   ("old_x_old_stock", "new_x_new_aug"),
    "vocab":  ("old_x_new_aug",   "new_x_new_aug"),
    "corpus": ("new_x_new_stock", "new_x_new_aug"),
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--hpoadj-jar", type=Path,
                        default=ANALYSIS_DIR.parents[2] / "hpoadj" / "hpoadj-cli" / "target" / "hpoadj-0.0.1.jar")
    parser.add_argument("--lirical-jar", type=Path,
                        default=ANALYSIS_DIR / "_tools" / "lirical-cli-2.4.1" / "lirical-cli-2.4.1.jar")
    parser.add_argument("--comparisons", nargs="+", choices=sorted(COMPARISONS), default=sorted(COMPARISONS),
                        help="which comparisons to (re)run; jobs shared between them run once")
    parser.add_argument("--cohorts", nargs="+", default=["APDS1", "NFKB1", "SOCS1"])
    parser.add_argument("--skip-adjust", action="store_true",
                        help="reuse existing hpoadj output where adjustment_summary.tsv exists")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def run(cmd, dry_run):
    print("  $", " ".join(str(c) for c in cmd))
    if not dry_run:
        subprocess.run([str(c) for c in cmd], check=True)


def kb_root(kb: str) -> Path:
    return WORK / "hpoa_adjusted" / KNOWLEDGE_BASES[kb]["outdir"]


def adjust_hpoa(kb: str, args):
    """Augment (if the knowledge base carries the cohort) and leave-one-publication-out."""
    spec = KNOWLEDGE_BASES[kb]
    outdir = kb_root(kb)
    if args.skip_adjust and (outdir / "adjustment_summary.tsv").exists():
        print(f"[adjust:{kb}] reusing {outdir}")
        return
    print(f"[adjust:{kb}] release={spec['release']} augment={spec['augment']}")
    hpoa = WORK / "data" / spec["release"] / "phenotype.hpoa"
    common = ["--hpo", WORK / "data" / spec["release"] / "hp.json", "-p", COHORT_DIR, "-o", outdir]
    if spec["augment"]:
        run(["java", "-jar", args.hpoadj_jar, "augment", "-a", hpoa,
             "--biocuration", "HPO:esid4hpo"] + common, args.dry_run)
        hpoa = outdir / "phenotype_augmented.hpoa"
    run(["java", "-jar", args.hpoadj_jar, "loo", "-a", hpoa] + common, args.dry_run)


def read_summary(kb: str):
    path = kb_root(kb) / "adjustment_summary.tsv"
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


def build_bundles(kb: str, adjusted, dry_run):
    """One LIRICAL data directory per publication: the stock bundle with the LOO phenotype.hpoa."""
    stock = WORK / "data" / KNOWLEDGE_BASES[kb]["release"]
    bundles = {}
    for pmid, hpoa in adjusted.items():
        bundle = kb_root(kb) / ("bundle_" + pmid.replace(":", "_"))
        bundles[pmid] = bundle
        if dry_run:
            print(f"[bundle:{kb}] would build {bundle}")
            continue
        if bundle.exists():
            shutil.rmtree(bundle)
        bundle.mkdir(parents=True)
        for entry in stock.iterdir():
            if entry.name != "phenotype.hpoa":
                (bundle / entry.name).symlink_to(entry.resolve())
        shutil.copy(hpoa, bundle / "phenotype.hpoa")
        print(f"[bundle:{kb}] {bundle.name} <- {hpoa.name}")
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


def benchmark(cohort: str, job: str, bundles, args):
    vocab, kb, merged_name = JOBS[job]
    staged = WORK / "lirical" / cohort / vocab
    packets = sorted(staged.glob("*.json"))
    if not packets:
        print(f"[lirical:{cohort}:{job}] no staged phenopackets in {staged}, skipping")
        return
    groups = defaultdict(list)
    for packet in packets:
        pmid = phenopacket_pmid(packet)
        groups[pmid if pmid in bundles else "stock"].append(packet)

    parts_dir = WORK / "lirical" / "parts"
    parts_dir.mkdir(parents=True, exist_ok=True)
    stock = WORK / "data" / KNOWLEDGE_BASES[kb]["release"]
    parts = []
    for label, members in sorted(groups.items(), key=lambda item: item[0]):
        data_dir = stock if label == "stock" else bundles[label]
        part = parts_dir / f"{cohort}.{job}.{label.replace(':', '_')}.csv"
        print(f"[lirical:{cohort}:{job}] {label}: {len(members)} phenopackets vs {data_dir.name}")
        run(["java", "-jar", args.lirical_jar, "benchmark",
             "-d", data_dir, "-o", part] + members, args.dry_run)
        parts.append(part)

    merged = WORK / "lirical" / merged_name.format(cohort=cohort)
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
    for jar in (args.hpoadj_jar, args.lirical_jar):
        if not args.dry_run and not jar.exists():
            sys.exit(f"jar not found: {jar}")

    jobs = list(dict.fromkeys(j for cmp in args.comparisons for j in COMPARISONS[cmp]))
    kbs = list(dict.fromkeys(JOBS[j][1] for j in jobs))
    print(f"[plan] comparisons={args.comparisons} -> jobs={jobs} -> knowledge bases={kbs}")

    for kb in kbs:
        adjust_hpoa(kb, args)

    bundles_by_kb, notes = {}, []
    for kb in kbs:
        adjusted, insufficient = read_summary(kb)
        bundles_by_kb[kb] = build_bundles(kb, adjusted, args.dry_run)
        if insufficient:
            seen = sorted({(pmid, disease) for _, pmid, disease in insufficient})
            notes.append(f"[note:{kb}] {len(insufficient)} phenopacket(s) have no phenotype annotation left "
                         "after removing their publication: "
                         + "; ".join(f"{pmid} is the sole source of {disease}" for pmid, disease in seen))

    for job in jobs:
        for cohort in args.cohorts:
            benchmark(cohort, job, bundles_by_kb[JOBS[job][1]], args)

    for note in notes:
        print(note)
    print("\nDone. Re-run the notebook from section 4 to rebuild rank tables and figures.")


if __name__ == "__main__":
    main()
