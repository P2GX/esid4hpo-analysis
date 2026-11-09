# esid4hpo-analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![ESID](https://img.shields.io/badge/ESID-European%20Society%20for%20Immunodeficiencies-blue)](https://esid.org/)
[![HPO](https://img.shields.io/badge/HPO-Human%20Phenotype%20Ontology-purple.svg)](https://hpo.jax.org/)
[![LIRICAL](https://img.shields.io/badge/LIRICAL-2.4.1-blue)](https://github.com/TheJacksonLaboratory/LIRICAL)

Analysis code for the ESID4HPO effort: quantifying what the ESID/HPO workshop added to the
Human Phenotype Ontology, by comparing cohorts of individuals annotated with two HPO releases.

| arm | release | how it is obtained |
|------|---------|--------------------|
| baseline | `v2024-08-13` | phenopackets are *aged* down to this vocabulary |
| workshop | `v2026-06-23` | phenopackets as curated |

Three cohorts are curated from primary literature with the post-workshop vocabulary:
- **APDS** (`OMIM:615513`, n = 38), 
- **NFKB1** (`OMIM:616576`, n = 52) and
- **SOCS1** (`OMIM:619375`, n = 28).

Two analyses are run over the same cohorts:

1. **Semantic content**: information content and specificity of the annotations, number of
   terms per individual, and the footprint of post-workshop terms.
2. **LIRICAL disease ranking**: rank of the causal disease per individual, before vs after.

## Repository layout

```
src/
  cohorts/<COHORT>/phenopackets/     curated phenopackets (primary data)
  cohorts/<COHORT>/*_individuals.json  cohort export from the curation tool
  analysis/
    esid4hpo_analysis.ipynb          main notebook: both analyses and all figures
    run_lirical_adjusted.py          HPOA augmentation + leave-one-publication-out + LIRICAL
    diagnose_profile_overlap.py      per-arm verbatim overlap between individuals and their LOO profiles
    profile_sizes.py                 phenotype annotations per disease profile and leave-one-out fold
    rerun_all.sh                     build hpoadj, run the above, print next steps
    zenodo/                          package the generated HPOA files for Zenodo, fetch them back
    ehr_ic.py                        EHR arm: IC per individual and source from the Southampton term counts
    figures/                         fig1_ontology_stats.py (Figure 1c), fig3_panels.py (run by the notebook), fig5_ehr.py (Figure 5, S2)
    _tools/                          LIRICAL distribution (not tracked)
    _work/                           outputs
      tables/                        manuscript and supplementary tables, named by supplement number
      intermediate/                  per-individual and per-term tables the summaries are built from
      figures/                       rendered figure panels (fig1C, fig3A, fig3B, fig4, fig5, figS1 to figS4)
      data/{old,new}/                version-matched LIRICAL data bundles (not tracked)
      hpoa_adjusted/{old,new_stock,new}/  generated HPOA files (not tracked, on Zenodo) + summary TSVs
      lirical/                       staged phenopackets and benchmark output (not tracked)
```

`fig3_panels.py` is executed from inside the notebook. `fig1_ontology_stats.py` computes the
term-level statistics of the immune branch between the two releases (Figure 1c, Supplementary
Table S2) and is run on its own from `src/analysis/figures/`; it needs the two `hp.json` files of
the data bundles. `ehr_ic.py` and `fig5_ehr.py` produce the EHR arm (Figure 5, Supplementary
Figure S2, Supplementary Table S5) from the Southampton term counts, which are not distributed.
The manuscript figure layouts themselves are assembled in PowerPoint and are not part of the
repository.

| Supplement | file in `src/analysis/_work/` |
| --- | --- |
| Table S2 (Figure 1c counts) | `tables/tableS2_immune_branch_counts.csv` |
| Table S3a to d (information content) | `tables/tableS3a_ic_summary.csv`, `tableS3b_ic_undefined.csv`, `tableS3c_term_novelty.csv`, `tableS3d_ic_structural.csv` |
| Table S4a to d (LIRICAL) | `tables/tableS4a_lirical_vocabulary.csv`, `tableS4b_lirical_annotations.csv`, `tableS4c_lirical_main.csv`, `tableS4d_lirical_main_loo.csv` |
| Table S5 (EHR) | `tables/tableS5_ehr_summary.csv` |
| Table S6 (excluded features omitted) | `tables/tableS6_ic_excluded_omitted.csv` |
| Table S7 (disease profile sizes) | `tables/tableS7_profile_sizes.csv` |
| Table S8a to c (separation) | `tables/tableS8a_separation_cohorts.csv`, `tableS8b_separation_sensitivity.csv`, `tableS8c_separation_publications.csv` |
| Figures 1c, 3a, 3b, 4, 5 | `figures/fig1C_immune_branch.*`, `fig3A_paired_ic.*`, `fig3B_post_workshop_fraction.*`, `fig4_lirical_decomposition.*`, `fig5_ehr.*` |
| Figures S1 to S4 | `figures/figS1_lirical_main.*`, `figS2_ehr_terms.*`, `figS3_ic_structural.*`, `figS4a_separation_mds.*`, `figS4b_separation_similarity.*` |

`figures/lirical_main_loo_sensitivity.*` is the figure counterpart of Table S4d and is not part
of the supplement. The hpoadj provenance tables are
`hpoa_adjusted/<knowledge base>/{augmentation,adjustment}_summary.tsv`.

## Why the HPOA was modified

The `v2026-06-23` HPOA release contains the workshop's new terms in `hp.json`, but not the
cohort curation: the annotation lines of the three benchmark diseases are identical to the
`v2024-08-13` release. Without ingesting the cohorts, the disease profiles do not differ
between arms and the comparison cannot detect any effect.

`run_lirical_adjusted.py` therefore calls [hpoadj](https://github.com/P2GX/hpoadj) in two
steps (`hpoadj augment`, `hpoadj loo`).

**Augmentation.** For each disease, one annotation per observed HPO term is written, with the
frequency pooled over all cohort publications using true-path (ancestor-aware) counting, the
denominator equal to the cohort size, and every contributing PMID in the reference field.
Pre-existing annotations are superseded where they derive from a publication the cohort
re-curates, and pooled with the cohort counts where they derive from an independent
publication. In the current run this raises the disease profiles to 149 (`OMIM:615513`), 285
(`OMIM:616576`) and 176 (`OMIM:619375`) distinct observed phenotype terms; the number of
annotation lines added, pooled and superseded per disease is written to
`_work/hpoa_adjusted/<arm>/augmentation_summary.tsv`. Only phenotype annotations (aspect `P`)
are superseded or pooled; onset and inheritance annotations are kept as released. A pre-existing
line without a frequency value is pooled as a 1/1 case report, which is how LIRICAL's annotation
loader (phenol, `HpoDiseaseLoaderDefault.parseFrequency`) reads such a line; negated,
percentage and frequency-term lines are kept beside the cohort line (none occur for the three
diseases). Excluded phenopacket features are not converted into `NOT` annotations, since absence
in a subset of individuals does not establish absence in the disease. The denominator is the
cohort size because phenopackets do not record which features were assessed.

**Leave-one-publication-out.** Phenopackets and disease annotations derive from the same
publications, so each individual must be scored against annotations that exclude its own source.
For every publication a copy of the HPOA is written in which that publication's contribution
is subtracted from the pooled counts (numerator and denominator reduced by its own individual
counts) and its identifier removed from the reference field; annotations resting solely on
that publication are dropped (`sole_support_removed_total` in `adjustment_summary.tsv`).
Only annotations pooled by hpoadj are decomposed; a curated annotation citing several
publications would be dropped and counted as `undecomposable_removed_total` (none occur for
the three diseases). The subtraction is exact: leaving every cohort publication out again
restores the released counts of a pooled annotation (hpoadj round-trip test). The biocuration
date of generated lines is fixed (`--biocuration-date`, default 2026-11-08) so that reruns
produce identical files.

## Reproducing the analysis

Prerequisites: Java 21, Maven, Python 3.12 (install the pinned packages with
`pip install -r src/analysis/requirements.txt`; the file is a `pip freeze` of the analysis
environment), the [hpoadj](https://github.com/P2GX/hpoadj) repository (v0.0.2) checked out
next to this one (the tool grew out of [hpotools PR #20](https://github.com/P2GX/hpotools/pull/20);
design discussion in issues [#17](https://github.com/P2GX/hpotools/issues/17),
[#18](https://github.com/P2GX/hpotools/issues/18) and
[#19](https://github.com/P2GX/hpotools/issues/19)), LIRICAL v2.4.1 unpacked in
`src/analysis/_tools/`, and the two data bundles prepared in `src/analysis/_work/data/{old,new}/`
(each with the release's `hp.json` and matching `phenotype.hpoa`, taken from the official HPO
release artifacts of `v2024-08-13` and `v2026-06-23`).

```bash
# 1. HPOA augmentation and LIRICAL benchmark (~80 min for all three comparisons;
#    `--comparisons main` for the headline pair only, ~40 min).
#    Needs the aged phenopackets staged by notebook section 3: on a fresh checkout,
#    first run notebook sections 0, 1 and 3 (section 2 needs the augmented HPOA
#    that this step produces, so skip it on the first pass).
cd src/analysis && bash rerun_all.sh

# 2. both notebooks, Restart & Run All: semantic analysis (needs the augmented
#    HPOA from step 1), LIRICAL rank tables, supplementary tables and figures,
#    then the separation analysis (Supplementary Figure S4, Table S8)
jupyter lab src/analysis/esid4hpo_analysis.ipynb
jupyter lab src/analysis/esid4hpo_separation.ipynb
```

`run_lirical_adjusted.py --dry-run` prints the plan without running anything. Two things
vary between arms: the vocabulary an individual is described in (pre-workshop terms from the
aged phenopackets, or post-workshop terms) and the knowledge base ranked against (the released
`phenotype.hpoa`, or the file augmented with the cohort annotations). The augmented file is
always leave-one-publication-out; the released v2024-08-13 file is used as released, which is
what a diagnostic laboratory had at the time, and additionally with the source publication
removed as a stricter sensitivity analysis. `--comparisons {main,main_loo,vocab,corpus}`
selects which paired comparisons to run (default: all four; jobs shared between them run once):

| comparison | pre arm | post arm | isolates |
| --- | --- | --- | --- |
| `main` | pre-workshop terms x released v2024-08-13 HPOA (`<cohort>.old_noloo.csv`) | post-workshop terms x augmented v2026-06-23 HPOA (`<cohort>.new.csv`) | the ESID4HPO output as a whole |
| `main_loo` | pre-workshop terms x released v2024-08-13 HPOA, source publication removed (`<cohort>.old.csv`) | as above | sensitivity analysis of `main` |
| `vocab` | pre-workshop terms x augmented v2026-06-23 HPOA (`<cohort>.old_terms_vs_new_hpoa.csv`) | as above | the new and re-parented terms alone |
| `corpus` | post-workshop terms x released v2026-06-23 HPOA (`<cohort>.new_terms_vs_stock_hpoa.csv`) | as above | the added disease annotations alone |

Knowledge bases are built under `_work/hpoa_adjusted/{old,new_stock,new}/`; `--skip-adjust`
reuses them. Notebook section 4 loads all five condition files, writes Supplementary Tables S4a
to S4d, Figure 4 (`_work/figures/fig4_lirical_decomposition.*`), Supplementary Figure S1
(`figS1_lirical_main.*`) and the `main_loo` figure (`lirical_main_loo_sensitivity.*`).
`profile_sizes.py` (run from the notebook) tabulates, per disease and left-out publication, the
phenotype annotations in every knowledge base an individual is ranked against (Table S7).
`diagnose_profile_overlap.py` reports, per arm, how many of each individual's observed terms
survive verbatim in the leave-one-publication-out profile the individual is scored against.

## What is tracked

Tracked: phenopackets, notebooks, scripts, the rendered figure panels (`_work/figures/`), the
manuscript and supplementary tables (`_work/tables/`), the intermediate tables they derive from
(`_work/intermediate/`, except the individual-level EHR tables) and the hpoadj provenance
summaries `_work/hpoa_adjusted/<knowledge base>/{augmentation,adjustment}_summary.tsv`.

Not tracked: the HPO/LIRICAL data bundles, the LIRICAL distribution, the generated HPOA
copies (~1 GB) and the full LIRICAL benchmark output (~300 MB). All of these are regenerated
by the commands above. The generated HPOA files (the augmented corpus and every
leave-one-publication-out copy, for all three knowledge bases) are deposited on Zenodo,
doi:10.5281/zenodo.24058794, and can be restored without rerunning hpoadj:

```bash
cd src/analysis && bash zenodo/fetch_hpoa_adjusted.sh
```

`zenodo/package_hpoa_for_zenodo.sh` builds the deposited archives and checksums from
`_work/hpoa_adjusted/`.
