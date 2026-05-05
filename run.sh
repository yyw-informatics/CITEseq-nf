#!/usr/bin/env bash
# Wrapper for `nextflow run` that redirects logs out of the project root.
# Usage: nextflow/run.sh [nextflow args...]
#   e.g. nextflow/run.sh nextflow/main.nf -profile sge -entry from_qc

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
LOG_DIR="${PROJECT_ROOT}/logs/nextflow"
mkdir -p "${LOG_DIR}"

export NXF_LOG_FILE="${LOG_DIR}/nextflow.log"

exec nextflow "$@"
