#!/usr/bin/env bash
# Packs the generated HPOA files under _work/hpoa_adjusted/ into one archive per knowledge
# base, plus README and checksums, ready for upload to Zenodo.
#
#   bash zenodo/package_hpoa_for_zenodo.sh            # from src/analysis
#
# Output: _work/zenodo/hpoa_adjusted_{new,new_stock,old}.tar.gz, README_hpoa_adjusted.md, SHA256SUMS
set -euo pipefail

ANALYSIS_DIR="$(cd "$(dirname "$0")/.." && pwd)"
SOURCE="${ANALYSIS_DIR}/_work/hpoa_adjusted"
OUT="${ANALYSIS_DIR}/_work/zenodo"

[[ -d "${SOURCE}" ]] || { echo "no generated files at ${SOURCE}; run run_lirical_adjusted.py --adjust-only first" >&2; exit 1; }

rm -rf "${OUT}"
mkdir -p "${OUT}"
export COPYFILE_DISABLE=1

for kb in new new_stock old; do
    archive="${OUT}/hpoa_adjusted_${kb}.tar.gz"
    tar -czf "${archive}" -C "${SOURCE}" --exclude='bundle_*' --exclude='.DS_Store' "${kb}"
    echo "$(du -h "${archive}" | cut -f1)  ${archive##*/}  ($(find "${SOURCE}/${kb}" -maxdepth 1 -name '*.hpoa' | wc -l | tr -d ' ') hpoa files)"
done

cp "${ANALYSIS_DIR}/zenodo/README_hpoa_adjusted.md" "${OUT}/"

cd "${OUT}"
if command -v sha256sum >/dev/null; then
    sha256sum hpoa_adjusted_*.tar.gz README_hpoa_adjusted.md > SHA256SUMS
else
    shasum -a 256 hpoa_adjusted_*.tar.gz README_hpoa_adjusted.md > SHA256SUMS
fi
cat SHA256SUMS
echo "upload the contents of ${OUT} to the Zenodo record"
