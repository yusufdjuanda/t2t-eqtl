#!/bin/bash
#SBATCH --job-name=t2t-asefeas
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=logs/06-compare/%x_%j.out
#SBATCH --error=logs/06-compare/%x_%j.err
# Step 6: feasibility of allele-specific expression at the changed genes (ase-feasibility.py).
#
#   scripts/submit.sh scripts/06-compare/ase-feasibility.sbatch <cohort>
set -euo pipefail
COHORT=${1:?cohort}
ROOT=${T2T_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && git rev-parse --show-toplevel)}}
source "$ROOT/scripts/lib/common.sh"
t2t_tmpdir "$SLURM_JOB_ID"
t2t_env bio   # STAR, salmon, samtools, ... from pixi.toml (published runs: cluster modules of the same versions)
echo "start $(date -u +%FT%TZ) $COHORT"
cd $ROOT
/usr/bin/time -v pixi run python scripts/06-compare/ase-feasibility.py --cohort $COHORT --threads $SLURM_CPUS_PER_TASK
echo "end $(date -u +%FT%TZ)"
