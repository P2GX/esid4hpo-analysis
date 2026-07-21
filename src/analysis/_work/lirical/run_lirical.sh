#!/usr/bin/env bash
set -euo pipefail

# NOTE: `lirical` is a shell ALIAS on this machine, and aliases do not work inside
# scripts. We therefore define a function that calls the jar by full path.
# (There is also no `module` command on macOS - that line was an HPC/Lmod artefact.)
LIRICAL_JAR="/Users/adamgraefe/Documents/git/esid4hpo-analysis/src/analysis/_tools/lirical-cli-2.4.1/lirical-cli-2.4.1.jar"
lirical() { java -jar "$LIRICAL_JAR" "$@"; }

WORKDIR="/Users/adamgraefe/Documents/git/esid4hpo-analysis/src/analysis/_work/lirical"
OLD_DATA="${WORKDIR}/../data/old"   # built with v2024-08-13
NEW_DATA="${WORKDIR}/../data/new"   # built with v2026-06-23

# Set FORCE=1 to recompute arms whose CSV already exists:
#   FORCE=1 bash run_lirical.sh
FORCE="${FORCE:-0}"

run_arm () {           # run_arm <cohort> <old|new> <data_dir>
  local cohort="$1" arm="$2" data="$3"
  local out="${WORKDIR}/${cohort}.${arm}.csv"
  local ppdir="${WORKDIR}/${cohort}/${arm}"

  if [[ ! -d "$ppdir" ]]; then
    echo "  [${cohort}/${arm}] no staged phenopackets - skipping"
    return 0
  fi
  if [[ -s "$out" && "$FORCE" != "1" ]]; then
    echo "  [${cohort}/${arm}] already done ($(basename "$out")) - skipping (FORCE=1 to redo)"
    return 0
  fi
  echo "  [${cohort}/${arm}] running $(ls "$ppdir"/*.json | wc -l | tr -d ' ') phenopackets ..."
  lirical benchmark -d "$data" -o "$out" "$ppdir"/*.json
}

for COHORT in SOCS1 APDS NFKB1; do
  [[ -d "${WORKDIR}/${COHORT}" ]] || continue
  echo "==== ${COHORT} ===="
  run_arm "$COHORT" old "$OLD_DATA"
  run_arm "$COHORT" new "$NEW_DATA"
done

echo
echo "Results in ${WORKDIR}:"
ls -1 "${WORKDIR}"/*.csv 2>/dev/null || echo "  (none yet)"
