#!/bin/bash
#SBATCH --job-name=t2t-gwaslook
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/06-compare/%x_%j.out
#SBATCH --error=logs/06-compare/%x_%j.err
# Step 6: GWAS Catalog lookup for the changed genes (gwas-lookup.py).
#
#   scripts/submit.sh scripts/06-compare/gwas-lookup.sbatch <cohort> [--background K]
set -euo pipefail
COHORT=${1:?cohort}
ROOT=${T2T_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && git rev-parse --show-toplevel)}}
source "$ROOT/scripts/lib/common.sh"
t2t_tmpdir "$SLURM_JOB_ID"
echo "start $(date -u +%FT%TZ) $COHORT"
cd $ROOT
/usr/bin/time -v pixi run python scripts/06-compare/gwas-lookup.py --cohort $COHORT "${@:2}"
echo "end $(date -u +%FT%TZ)"
