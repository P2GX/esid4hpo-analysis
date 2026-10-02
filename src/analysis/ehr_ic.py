#!/usr/bin/env python3
"""Information content of HPO terms recovered from Southampton EHR text, pre- vs post-workshop.

Input: a patient-level table of concept-recognition hits (FastHPOCR), one row per
(patient, source, HPO term, HPO release), columns anon_id, source, HPO_ID, HPO_term,
count, HPO_version. The file is NOT part of the repository; by default it is read from
_work/data/ehr/hpo_term_counts.xlsx (gitignored).

Design. The same text was annotated once with HPO v2024-08-13 and once with v2026-06-23,
so no ageing step is needed. Terms are scored with cohort-based Resnik information content,
the same formula as the curated-cohort analysis with patients in place of diseases:

    IC(t) = -ln(freq(t) / N) / ln N

where freq(t) is the number of patients in the unit whose POST-workshop term set contains
t or any descendant of t (ancestor closure on the v2026-06-23 graph, restricted to the
Phenotypic abnormality subtree), and N is the number of patients in the unit. One frequency
table per unit scores BOTH arms, so a term present in both releases scores identically and
the paired difference reflects only the terms recovered.

Units: each data source on its own (icd10, clinic_letters, histology, radiology) and
"combined" = histology and radiology term sets unioned per patient.

Per patient and arm: number of distinct immune-branch terms, mean term IC, summed term IC.
Statistics: two-sided Wilcoxon signed-rank, post vs pre, on per-patient mean IC (patients
with at least one scorable term in BOTH arms) and on per-patient sum IC (all N patients,
empty set = 0); uncorrected, as in Alex's definitions.

Outputs (in _work/):
    ehr_ic_summary.csv         one row per unit (the manuscript table)       -> tracked
    ehr_ic_patient_level.csv   per patient, unit and arm                     -> gitignored
    ehr_ic_term_level.csv      IC per term and unit, with patient counts     -> tracked

Usage, from src/analysis:
    python ehr_ic.py
    python ehr_ic.py --counts /path/to/hpo_term_counts.xlsx --hpo _work/data/new/hp.json
"""

import argparse
import csv
import json
import math
import pathlib
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

try:
    from scipy.stats import wilcoxon as _wilcoxon
except Exception:  # pragma: no cover
    _wilcoxon = None

ANALYSIS = pathlib.Path(__file__).resolve().parent
WORK = ANALYSIS / "_work"
OLD_TAG, NEW_TAG = "2024-08-13", "2026-06-23"
ROOT = "HP:0000118"            # Phenotypic abnormality
IMMUNE = "HP:0002715"          # Abnormality of the immune system
UNITS = {                      # unit -> sources unioned per patient
    "icd10": ["icd10"],
    "clinic_letters": ["clinic_letters"],
    "histology": ["histology"],
    "imaging": ["radiology"],
    "combined": ["histology", "radiology"],
}


# ----------------------------------------------------------------------------- ontology
def load_hpo(path):
    """Return (parents, alt_id -> primary id, deprecated set, labels) from an obographs hp.json."""
    g = json.load(open(path))["graphs"][0]
    parents, alt, deprecated, label = defaultdict(list), {}, set(), {}
    for n in g["nodes"]:
        cid = n["id"].split("/")[-1].replace("_", ":")
        if not cid.startswith("HP:"):
            continue
        label[cid] = n.get("lbl", "")
        meta = n.get("meta", {})
        if meta.get("deprecated"):
            deprecated.add(cid)
        for bpv in meta.get("basicPropertyValues", []):
            if bpv.get("pred", "").endswith("hasAlternativeId"):
                alt[bpv["val"]] = cid
            if bpv.get("pred", "").endswith("IAO_0100001"):       # term replaced by
                alt.setdefault(cid, bpv["val"].split("/")[-1].replace("_", ":"))
    for e in g["edges"]:
        if e["pred"] == "is_a":
            parents[e["sub"].split("/")[-1].replace("_", ":")].append(
                e["obj"].split("/")[-1].replace("_", ":"))
    return parents, alt, deprecated, label


def make_ancestor_fn(parents, root=ROOT):
    cache = {}

    def ancestors(curie):
        if curie in cache:
            return cache[curie]
        seen, stack = {curie}, list(parents.get(curie, []))
        while stack:
            x = stack.pop()
            if x not in seen:
                seen.add(x)
                stack.extend(parents.get(x, []))
        out = seen if root in seen else set()
        cache[curie] = out
        return out

    return ancestors


# ----------------------------------------------------------------------------- IC
def ic_table(term_sets, ancestors):
    """Cohort-based Resnik IC: freq(t) over patients with t or a descendant; normalised by ln N."""
    n = len(term_sets)
    freq = defaultdict(int)
    for terms in term_sets:
        closure = set()
        for t in terms:
            closure |= ancestors(t)
        for t in closure:
            freq[t] += 1
    mx = math.log(n) if n > 1 else 1.0
    return {t: -math.log(min(f / n, 1.0)) / mx for t, f in freq.items()}, dict(freq), n


def wilcoxon_p(a, b):
    """Two-sided Wilcoxon signed-rank p for paired samples a (pre) and b (post)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = b - a
    d = d[d != 0]
    if len(d) == 0:
        return float("nan")
    if _wilcoxon is not None:
        return float(_wilcoxon(b, a).pvalue)
    # normal approximation with tie correction (fallback when scipy is missing)
    r = pd.Series(np.abs(d)).rank().values
    w = r[d > 0].sum()
    n = len(d)
    mu, var = n * (n + 1) / 4, n * (n + 1) * (2 * n + 1) / 24
    _, counts = np.unique(r, return_counts=True)
    var -= (counts ** 3 - counts).sum() / 48
    z = (w - mu) / math.sqrt(var)
    return float(math.erfc(abs(z) / math.sqrt(2)))


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", type=pathlib.Path, default=WORK / "data" / "ehr" / "hpo_term_counts.xlsx")
    ap.add_argument("--hpo", type=pathlib.Path, default=WORK / "data" / "new" / "hp.json")
    ap.add_argument("--out", type=pathlib.Path, default=WORK)
    args = ap.parse_args()

    parents, alt, deprecated, label = load_hpo(args.hpo)
    ancestors = make_ancestor_fn(parents)

    df = pd.read_excel(args.counts) if args.counts.suffix == ".xlsx" else pd.read_csv(args.counts)
    df["HPO_version"] = pd.to_datetime(df["HPO_version"]).dt.strftime("%Y-%m-%d")
    df["term"] = df["HPO_ID"].map(lambda t: alt.get(t, t))          # obsolete ids -> primary id
    df = df[df["term"].isin(parents.keys() | set(label))]
    unmapped = sorted(set(df.loc[df["term"].map(lambda t: not ancestors(t)), "term"]))
    if unmapped:
        print(f"[note] {len(unmapped)} term(s) outside the Phenotypic abnormality subtree or "
              f"deprecated without replacement are unscorable: {', '.join(unmapped[:10])}", file=sys.stderr)

    patient_rows, term_rows, summary = [], [], []
    for unit, sources in UNITS.items():
        sub = df[df["source"].isin(sources)]
        sets = {arm: sub[sub.HPO_version == tag].groupby("anon_id")["term"].agg(set).to_dict()
                for arm, tag in (("pre", OLD_TAG), ("post", NEW_TAG))}
        patients = sorted(set(sets["pre"]) | set(sets["post"]))
        n_unit = len(patients)
        # ONE frequency table per unit, from the post-workshop term sets, scores both arms
        table, freq, _ = ic_table([sets["post"].get(p, set()) for p in patients], ancestors)

        per = {}
        for p in patients:
            row = dict(unit=unit, anon_id=p)
            for arm in ("pre", "post"):
                terms = sets[arm].get(p, set())
                ics = [table[t] for t in terms if t in table]
                row[f"n_terms_{arm}"] = len(terms)
                row[f"n_unscored_{arm}"] = len(terms) - len(ics)
                row[f"mean_ic_{arm}"] = float(np.mean(ics)) if ics else np.nan
                row[f"sum_ic_{arm}"] = float(np.sum(ics))
            per[p] = row
            patient_rows.append(row)
        pdf = pd.DataFrame(list(per.values()))
        paired = pdf.dropna(subset=["mean_ic_pre", "mean_ic_post"])

        summary.append({
            "unit": unit, "N": n_unit, "paired_n": len(paired),
            "terms_per_patient_pre": pdf.n_terms_pre.mean(), "terms_per_patient_post": pdf.n_terms_post.mean(),
            "mean_ic_pre": paired.mean_ic_pre.mean(), "mean_ic_post": paired.mean_ic_post.mean(),
            "mean_ic_delta": (paired.mean_ic_post - paired.mean_ic_pre).mean(),
            "mean_ic_improved": int((paired.mean_ic_post > paired.mean_ic_pre).sum()),
            "mean_ic_worse": int((paired.mean_ic_post < paired.mean_ic_pre).sum()),
            "sum_ic_pre": pdf.sum_ic_pre.mean(), "sum_ic_post": pdf.sum_ic_post.mean(),
            "sum_ic_delta": (pdf.sum_ic_post - pdf.sum_ic_pre).mean(),
            "sum_ic_improved": int((pdf.sum_ic_post > pdf.sum_ic_pre).sum()),
            "sum_ic_worse": int((pdf.sum_ic_post < pdf.sum_ic_pre).sum()),
            "p_mean_ic": wilcoxon_p(paired.mean_ic_pre, paired.mean_ic_post),
            "p_sum_ic": wilcoxon_p(pdf.sum_ic_pre, pdf.sum_ic_post),
            "distinct_terms_pre": len(set().union(*sets["pre"].values())) if sets["pre"] else 0,
            "distinct_terms_post": len(set().union(*sets["post"].values())) if sets["post"] else 0,
            "unscored_hits_pre": int(pdf.n_unscored_pre.sum()), "unscored_hits_post": int(pdf.n_unscored_post.sum()),
        })

        # term level: IC and how many patients carry the term directly, per arm
        direct = {arm: defaultdict(int) for arm in ("pre", "post")}
        for arm in ("pre", "post"):
            for terms in sets[arm].values():
                for t in terms:
                    direct[arm][t] += 1
        for t in sorted(set(direct["pre"]) | set(direct["post"])):
            term_rows.append(dict(unit=unit, term=t, label=label.get(t, ""),
                                  ic=table.get(t, np.nan), freq_closure=freq.get(t, 0),
                                  patients_pre=direct["pre"].get(t, 0), patients_post=direct["post"].get(t, 0),
                                  immune_branch=IMMUNE in ancestors(t)))

    args.out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(summary).to_csv(args.out / "ehr_ic_summary.csv", index=False)
    pd.DataFrame(patient_rows).to_csv(args.out / "ehr_ic_patient_level.csv", index=False)
    pd.DataFrame(term_rows).to_csv(args.out / "ehr_ic_term_level.csv", index=False)

    s = pd.DataFrame(summary)
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print(s.round(3).to_string(index=False))
    print(f"\nwrote {args.out / 'ehr_ic_summary.csv'}, ehr_ic_patient_level.csv (do not commit), ehr_ic_term_level.csv")


if __name__ == "__main__":
    main()
