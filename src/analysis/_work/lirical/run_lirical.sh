#!/usr/bin/env bash
set -euo pipefail
module load lirical/2.2.1   # or: LIRICAL=path/to/lirical
WORKDIR="/Users/adamgraefe/Documents/git/esid4hpo-analysis/src/analysis/_work/lirical"
OLD_DATA="${WORKDIR}/../data/old"   # built with v2024-08-13
NEW_DATA="${WORKDIR}/../data/new"   # built with v2026-06-23

# ---- SOCS1 ----
lirical benchmark -d "${OLD_DATA}" -o "${WORKDIR}/SOCS1.old.csv" "${WORKDIR}/SOCS1/old"/*.json
lirical benchmark -d "${NEW_DATA}" -o "${WORKDIR}/SOCS1.new.csv" "${WORKDIR}/SOCS1/new"/*.json

# ---- APDS ----
lirical benchmark -d "${OLD_DATA}" -o "${WORKDIR}/APDS.old.csv" "${WORKDIR}/APDS/old"/*.json
lirical benchmark -d "${NEW_DATA}" -o "${WORKDIR}/APDS.new.csv" "${WORKDIR}/APDS/new"/*.json

# ---- NFKB1 ----
lirical benchmark -d "${OLD_DATA}" -o "${WORKDIR}/NFKB1.old.csv" "${WORKDIR}/NFKB1/old"/*.json
lirical benchmark -d "${NEW_DATA}" -o "${WORKDIR}/NFKB1.new.csv" "${WORKDIR}/NFKB1/new"/*.json
