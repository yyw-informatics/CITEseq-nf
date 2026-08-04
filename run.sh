#!/usr/bin/env bash
# Wrapper for `nextflow run` that keeps .nextflow.log* out of the working tree.
# Usage: ./run.sh [nextflow args...]
#   e.g. ./run.sh run ./main.nf -profile sge -entry from_qc
# Override the log location with CITESEQ_LOG_DIR.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
LOG_DIR="${CITESEQ_LOG_DIR:-${REPO_ROOT}/logs/nextflow}"
mkdir -p "${LOG_DIR}"

export NXF_LOG_FILE="${LOG_DIR}/nextflow.log"

exec nextflow "$@"
