#!/bin/bash
#SBATCH --job-name=t2t-phenocov
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=/storage/projects/yusuf.maulana/t2t-eqtl/logs/05-eqtl/%x_%j.out
#SBATCH --error=/storage/projects/yusuf.maulana/t2t-eqtl/logs/05-eqtl/%x_%j.err
# Step 5: phenotype BED and covariates for one cohort and quantification method, on all three references.
#
#   sbatch scripts/05-eqtl/prepare-phenotypes-covariates.sbatch <cohort> [method=star]
#
# Runs prepare-phenotypes.py then make-covariates.py per reference (genotype PCs are reused from the cohort cache).
set -euo pipefail
COHORT=${1:?cohort}; METHOD=${2:-star}
ROOT=/storage/home/yusuf.maulana/research/t2t-eqtl
export TMPDIR=/scratch/$USER/t2t-eqtl/tmp/$SLURM_JOB_ID; mkdir -p $TMPDIR; trap 'rm -rf $TMPDIR' EXIT
echo "start $(date -u +%FT%TZ) $COHORT $METHOD"
cd $ROOT
for r in chm13v2.0 grch38 grch38-ffx; do
  /usr/bin/time -v pixi run python scripts/05-eqtl/prepare-phenotypes.py --cohort $COHORT --ref $r --method $METHOD
  pixi run python scripts/05-eqtl/make-covariates.py --cohort $COHORT --ref $r --method $METHOD
done
echo "end $(date -u +%FT%TZ)"
