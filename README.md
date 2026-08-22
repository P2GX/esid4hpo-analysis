# esid4hpo-analysis

Analysis code for the ESID4HPO effort: quantifying what the ESID/HPO workshop added to the
Human Phenotype Ontology, by comparing patient cohorts annotated with two HPO releases.

| arm | release | how it is obtained |
|------|---------|--------------------|
| baseline | `v2024-08-13` | phenopackets are *aged* down to this vocabulary |
| workshop | `v2026-06-23` | phenopackets as curated |

Three cohorts are curated from primary literature with the post-workshop vocabulary:
- **APDS** (`OMIM:615513`, n = 38), 
- **NFKB1** (`OMIM:616576`, n = 24) and
- **SOCS1** (`OMIM:619375`, n = 28).

Two analyses are run over the same cohorts:

1. **Semantic content** — information content and specificity of the annotations, number of
   terms per patient, and the footprint of post-workshop terms.
2. **LIRICAL disease ranking** — rank of the causal disease per patient, before vs after.

## Repository layout

```
src/
  cohorts/<COHORT>/phenopackets/     curated phenopackets (primary data)
  cohorts/<COHORT>/*_individuals.json  cohort export from the curation tool
  analysis/
    esid4hpo_analysis.ipynb          main notebook: both analyses and all figures
    run_lirical_adjusted.py          HPOA augmentation + leave-one-publication-out + LIRICAL
    diagnose_profile_overlap.py      per-arm verbatim overlap between patients and their LOO profiles
    rerun_all.sh                     build hpotools, run the above, print next steps
    figures/                         standalone figure scripts
    _tools/                          LIRICAL distribution (not tracked)
    _work/                           scratch and outputs
      data/{old,new}/                version-matched LIRICAL data bundles (not tracked)
      hpoa_adjusted/{old,new}/       generated HPOA files (not tracked) + summary TSVs
      lirical/                       benchmark output (not tracked) + staged phenopackets
      figures/                       rendered manuscript figures
      *.csv                          result tables
```

## Why the HPOA was modified

The `v2026-06-23` HPOA release contains the workshop's new terms in `hp.json`, but not the
cohort curation: the annotation lines of the three benchmark diseases are identical to the
`v2024-08-13` release. Without ingesting the cohorts, the disease profiles do not differ
between arms and the comparison cannot detect any effect.

`run_lirical_adjusted.py` therefore calls the `hpoadjust` command of
[hpotools](https://github.com/P2GX/hpotools) in two steps.

**Augmentation.** For each disease, one annotation per observed HPO term is written, with the
frequency pooled over all cohort publications using true-path (ancestor-aware) counting, the
denominator equal to the cohort size, and every contributing PMID in the reference field.
Pre-existing annotations are superseded where they derive from a publication the cohort
re-curates, and pooled with the cohort counts where they derive from an independent
publication. In the current run this adds 149 (`OMIM:615513`), 103 (`OMIM:616576`) and 180
(`OMIM:619375`) pooled annotations; the exact per-disease counts of every run are written to
`_work/hpoa_adjusted/<arm>/augmentation_summary.tsv`.

**Leave-one-publication-out.** Phenopackets and disease annotations derive from the same
publications, so each patient must be scored against annotations that exclude its own source.
For every publication a copy of the HPOA is written in which that publication's contribution
is subtracted from the pooled counts (numerator and denominator reduced by its own patient
counts) and its identifier removed from the reference field; annotations resting solely on
that publication are dropped. The subtraction is exact: reversing the cohort contribution of
a pooled annotation restores the original HPOA counts.

## Reproducing the analysis

Prerequisites: Java 21, Maven, Python 3.11+ (install the pinned packages with
`pip install -r src/analysis/requirements.txt`; the file is generated with `pip freeze`
from the analysis environment), the [hpotools](https://github.com/P2GX/hpotools) repository
checked out next to this one on branch `hpoa-adjuster-module` at commit `c7101e1`
([PR #20](https://github.com/P2GX/hpotools/pull/20), which adds the `hpoadjust` command;
design discussion in issues [#17](https://github.com/P2GX/hpotools/issues/17),
[#18](https://github.com/P2GX/hpotools/issues/18) and
[#19](https://github.com/P2GX/hpotools/issues/19)), LIRICAL v2.4.1 unpacked in
`src/analysis/_tools/`, and the two data bundles prepared in `src/analysis/_work/data/{old,new}/`
(each with the release's `hp.json` and matching `phenotype.hpoa`, taken from the official HPO
release artifacts of `v2024-08-13` and `v2026-06-23`).

```bash
# 1. sections 1-3 of the notebook: aged phenopackets, semantic analysis, LIRICAL inputs
jupyter lab src/analysis/esid4hpo_analysis.ipynb

# 2. HPOA augmentation and LIRICAL benchmark (~40 min, 36 LIRICAL runs)
cd src/analysis && bash rerun_all.sh

# 3. sections 4+ of the notebook: rank tables, statistics, figures
```

`run_lirical_adjusted.py --dry-run` prints the plan without running anything.
`--augment-arms {new,both,none}` selects which arm receives the cohort annotations;
`new` (default) compares the published baseline against the workshop release plus curation,
`both` holds the curation constant and isolates the effect of the vocabulary alone.
Non-default modes write to suffixed paths (`_work/hpoa_adjusted.<mode>/`,
`_work/lirical/<cohort>.<arm>.<mode>.csv`), so a decomposition run cannot overwrite the
headline results. `diagnose_profile_overlap.py` reports, per arm, how many of each patient's
observed terms survive verbatim in the leave-one-publication-out profile the patient is
scored against.

## What is tracked

Tracked: phenopackets, notebook, scripts, rendered figures, and the small result tables
(`_work/lirical_ranks_long.csv`, `_work/lirical_ranks_long.both.csv`, `_work/semantic_*.csv`,
`_work/fig3_ic_stats.csv`) plus the provenance summaries
`_work/hpoa_adjusted*/<arm>/{augmentation,adjustment}_summary.tsv` for every augmentation mode.

Not tracked: the HPO/LIRICAL data bundles, the LIRICAL distribution, the generated HPOA
copies (~1.3 GB) and the full LIRICAL benchmark output (~300 MB). All of these are
regenerated by the commands above; the generated HPOA files will be deposited alongside the
manuscript.
