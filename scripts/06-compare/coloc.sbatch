#!/bin/bash
#SBATCH --job-name=t2t-coloc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=/storage/projects/yusuf.maulana/t2t-eqtl/logs/06-compare/%x_%j.out
#SBATCH --error=/storage/projects/yusuf.maulana/t2t-eqtl/logs/06-compare/%x_%j.err
# Step 6 (H6): colocalisation of the changed genes with GWAS on each reference (coloc.py).
#
#   sbatch scripts/06-compare/coloc.sbatch <cohort> [--traits GCST...] [--out-suffix S]
set -euo pipefail
COHORT=${1:?cohort}
ROOT=/storage/home/yusuf.maulana/research/t2t-eqtl
export TMPDIR=/scratch/$USER/t2t-eqtl/tmp/$SLURM_JOB_ID; mkdir -p $TMPDIR; trap 'rm -rf $TMPDIR' EXIT
echo "start $(date -u +%FT%TZ) $COHORT"
cd $ROOT
/usr/bin/time -v pixi run python scripts/06-compare/coloc.py --cohort $COHORT "${@:2}"
echo "end $(date -u +%FT%TZ)"
