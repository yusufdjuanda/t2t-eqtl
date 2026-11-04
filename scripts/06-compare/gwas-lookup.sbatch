#!/bin/bash
#SBATCH --job-name=t2t-gwaslook
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/storage/projects/yusuf.maulana/t2t-eqtl/logs/06-compare/%x_%j.out
#SBATCH --error=/storage/projects/yusuf.maulana/t2t-eqtl/logs/06-compare/%x_%j.err
# Step 6: GWAS Catalog lookup for the changed genes (gwas-lookup.py).
#
#   sbatch scripts/06-compare/gwas-lookup.sbatch <cohort> [--background K]
set -euo pipefail
COHORT=${1:?cohort}
ROOT=/storage/home/yusuf.maulana/research/t2t-eqtl
export TMPDIR=/scratch/$USER/t2t-eqtl/tmp/$SLURM_JOB_ID; mkdir -p $TMPDIR; trap 'rm -rf $TMPDIR' EXIT
echo "start $(date -u +%FT%TZ) $COHORT"
cd $ROOT
/usr/bin/time -v pixi run python scripts/06-compare/gwas-lookup.py --cohort $COHORT "${@:2}"
echo "end $(date -u +%FT%TZ)"
