#!/bin/bash
#SBATCH --job-name=t2t-merge
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=3
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/storage/projects/yusuf.maulana/t2t-eqtl/logs/04-expression/%x_%j.out
#SBATCH --error=/storage/projects/yusuf.maulana/t2t-eqtl/logs/04-expression/%x_%j.err
# Stage 04 (expression): gene x individual matrices for one cohort (merge-expression.py), one process per reference.
#
#   sbatch scripts/04-expression/merge-expression.sbatch <cohort>
set -euo pipefail
COHORT=${1:?cohort}
ROOT=/storage/home/yusuf.maulana/research/t2t-eqtl
export TMPDIR=/scratch/$USER/t2t-eqtl/tmp/$SLURM_JOB_ID; mkdir -p $TMPDIR; trap 'rm -rf $TMPDIR' EXIT
echo "start $(date -u +%FT%TZ) $COHORT"
cd $ROOT
/usr/bin/time -v pixi run python scripts/04-expression/merge-expression.py --cohort $COHORT --jobs 3
echo "end $(date -u +%FT%TZ)"
