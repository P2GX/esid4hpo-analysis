# ESID4HPO: adjusted HPO annotation corpus files

Annotation corpus files (`phenotype.hpoa` format) generated for the ESID4HPO evaluation of
expert HPO curation in inborn errors of immunity. They are the knowledge bases against which
LIRICAL ranked the 118 published individuals with SOCS1 haploinsufficiency (OMIM:619375),
activated PI3Kδ syndrome 1 (OMIM:615513) and NFKB1 haploinsufficiency (OMIM:616576), and the
frequency tables of the information-content analysis.

Generated with hpoadj v0.0.2 (https://github.com/P2GX/hpoadj, doi:10.5281/zenodo.24054540)
by `run_lirical_adjusted.py` in https://github.com/P2GX/esid4hpo-analysis, from the cohort
Phenopackets in that repository and the HPO release artefacts listed below. The biocuration
date of generated lines is fixed to 2026-11-08, so a rerun reproduces these files byte for byte.

## Archives

| archive | release | contents |
| --- | --- | --- |
| `hpoa_adjusted_new.tar.gz` | v2026-06-23 | `phenotype_augmented.hpoa`: released corpus augmented with the three cohorts; `phenotype_PMID_<n>.hpoa`: one leave-one-publication-out copy per source publication (20) |
| `hpoa_adjusted_new_stock.tar.gz` | v2026-06-23 | released corpus, unaugmented, leave-one-publication-out for the three publications it cites (`PMID:24136356`, `PMID:24165795`, `PMID:33087723`) |
| `hpoa_adjusted_old.tar.gz` | v2024-08-13 | released corpus of the pre-workshop release, leave-one-publication-out for the same three publications |

Each archive also holds the hpoadj provenance tables: `augmentation_summary.tsv` (per disease:
cohort size, PMIDs, annotation lines added, pooled, superseded, retained beside the cohort line)
and `adjustment_summary.tsv` (per phenopacket: lines removed and subtracted, sole-support and
undecomposable lines, phenotype annotations remaining in the disease profile). `SHA256SUMS`
lists the checksums of the archives.

## Source data

| input | version |
| --- | --- |
| `hp.json`, `phenotype.hpoa` (post-workshop) | HPO release v2026-06-23 |
| `hp.json`, `phenotype.hpoa` (pre-workshop) | HPO release v2024-08-13 |
| cohort Phenopackets (GA4GH Phenopacket Schema v2) | esid4hpo-analysis, `src/cohorts/` |

## How the files were built

For every HPO term directly observed in a cohort, one annotation line was written whose
frequency is the number of individuals carrying the term or any of its descendants, divided by
the cohort size, with every cohort PMID in the reference field. Pre-existing phenotype
annotations were superseded where their publication was re-curated by the cohort and pooled
with the cohort counts where it was not. For each source publication, a copy was derived in
which that publication's counts are subtracted from numerators and denominators, its PMID
removed from the reference field, and lines supported by it alone dropped. Onset and
inheritance annotations and all other diseases are unchanged from the release.

## Licence

Derived from the HPO annotation files, which are distributed under the HPO licence
(https://hpo.jax.org/license). The derived files are made available under the same terms.
