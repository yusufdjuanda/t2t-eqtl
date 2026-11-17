#!/bin/bash
#SBATCH --job-name=t2t-h1qc
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=logs/06-compare/%x_%j.out
#SBATCH --error=logs/06-compare/%x_%j.err
# Step 6: quality check of the H1 lead variants (h1-variant-qc.py) for one cohort and quantification method.
#
#   scripts/submit.sh scripts/06-compare/h1-variant-qc.sbatch [cohort=mage] [method=star]
set -euo pipefail
COHORT=${1:-mage}; METHOD=${2:-star}
ROOT=${T2T_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && git rev-parse --show-toplevel)}}
source "$ROOT/scripts/lib/common.sh"
t2t_tmpdir "$SLURM_JOB_ID"
t2t_env bio   # STAR, salmon, samtools, ... from pixi.toml (published runs: cluster modules of the same versions)
echo "start $(date -u +%FT%TZ) $COHORT $METHOD"
cd $ROOT
/usr/bin/time -v pixi run python scripts/06-compare/h1-variant-qc.py --cohort $COHORT --method $METHOD
echo "end $(date -u +%FT%TZ)"
