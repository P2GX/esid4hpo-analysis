#!/usr/bin/env bash
# Downloads the generated HPOA files from Zenodo and unpacks them into _work/hpoa_adjusted/,
# the location run_lirical_adjusted.py and the notebook read from.
#
#   bash zenodo/fetch_hpoa_adjusted.sh            # from src/analysis
#   ZENODO_RECORD=1234567 bash zenodo/fetch_hpoa_adjusted.sh
set -euo pipefail

ZENODO_RECORD="${ZENODO_RECORD:-24058794}"
BASE="https://zenodo.org/records/${ZENODO_RECORD}/files"

ANALYSIS_DIR="$(cd "$(dirname "$0")/.." && pwd)"
TARGET="${ANALYSIS_DIR}/_work/hpoa_adjusted"
DOWNLOAD="${ANALYSIS_DIR}/_work/zenodo_download"

mkdir -p "${DOWNLOAD}" "${TARGET}"
cd "${DOWNLOAD}"

for name in SHA256SUMS README_hpoa_adjusted.md hpoa_adjusted_new.tar.gz hpoa_adjusted_new_stock.tar.gz hpoa_adjusted_old.tar.gz; do
    echo "fetching ${name}"
    curl -fsSL -o "${name}" "${BASE}/${name}?download=1"
done

if command -v sha256sum >/dev/null; then
    sha256sum -c SHA256SUMS
else
    shasum -a 256 -c SHA256SUMS
fi

for kb in new new_stock old; do
    tar -xzf "hpoa_adjusted_${kb}.tar.gz" -C "${TARGET}"
done
echo "unpacked into ${TARGET}"
