#!/usr/bin/env bash
set -euo pipefail

ANALYSIS_DIR="$(cd "$(dirname "$0")" && pwd)"
HPOTOOLS_DIR="${ANALYSIS_DIR}/../../../hpotools"

echo "== 1/3 build hpotools =="
(cd "${HPOTOOLS_DIR}" && mvn -q test && mvn -q -DskipTests package)

echo "== 2/3 augment the HPOA with the cohort annotations and run LIRICAL =="
echo "   36 LIRICAL runs (18 publications x 2 arms), roughly 40 minutes"
(cd "${ANALYSIS_DIR}" && python run_lirical_adjusted.py "$@")

echo "== 3/3 done =="
echo "Now open esid4hpo_analysis.ipynb and Restart & Run All."
echo "Figures land in _work/figures/, tables in _work/*.csv,"
echo "per-disease augmentation counts in _work/hpoa_adjusted/<arm>/augmentation_summary.tsv"
