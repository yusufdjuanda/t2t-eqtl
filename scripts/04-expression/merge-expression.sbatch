#!/bin/bash
#SBATCH --job-name=t2t-merge
#SBATCH --cpus-per-task=3
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=logs/04-expression/%x_%j.out
#SBATCH --error=logs/04-expression/%x_%j.err
# Stage 04 (expression): gene x individual matrices for one cohort (merge-expression.py), one process per reference.
#
#   scripts/submit.sh scripts/04-expression/merge-expression.sbatch <cohort>
set -euo pipefail
COHORT=${1:?cohort}
ROOT=${T2T_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && git rev-parse --show-toplevel)}}
source "$ROOT/scripts/lib/common.sh"
t2t_tmpdir "$SLURM_JOB_ID"
echo "start $(date -u +%FT%TZ) $COHORT"
cd $ROOT
/usr/bin/time -v pixi run python scripts/04-expression/merge-expression.py --cohort $COHORT --jobs 3
echo "end $(date -u +%FT%TZ)"
