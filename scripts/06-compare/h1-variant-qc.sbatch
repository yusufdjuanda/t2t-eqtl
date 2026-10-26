#!/bin/bash
#SBATCH --job-name=t2t-h1qc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=/storage/projects/yusuf.maulana/t2t-eqtl/logs/06-compare/%x_%j.out
#SBATCH --error=/storage/projects/yusuf.maulana/t2t-eqtl/logs/06-compare/%x_%j.err
# Step 6: quality check of the H1 lead variants (h1-variant-qc.py) for one cohort and quantification method.
#
#   sbatch scripts/06-compare/h1-variant-qc.sbatch [cohort=mage] [method=star]
set -euo pipefail
COHORT=${1:-mage}; METHOD=${2:-star}
ROOT=/storage/home/yusuf.maulana/research/t2t-eqtl
export TMPDIR=/scratch/$USER/t2t-eqtl/tmp/$SLURM_JOB_ID; mkdir -p $TMPDIR; trap 'rm -rf $TMPDIR' EXIT
module load bcftools/1.23.1
echo "start $(date -u +%FT%TZ) $COHORT $METHOD"
cd $ROOT
/usr/bin/time -v pixi run python scripts/06-compare/h1-variant-qc.py --cohort $COHORT --method $METHOD
echo "end $(date -u +%FT%TZ)"
