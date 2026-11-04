#!/bin/bash
#SBATCH --job-name=t2t-collect
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/storage/projects/yusuf.maulana/t2t-eqtl/logs/05-eqtl/%x_%j.out
#SBATCH --error=/storage/projects/yusuf.maulana/t2t-eqtl/logs/05-eqtl/%x_%j.err
# Step 5: per-chromosome tensorQTL results -> cis-eqtl.tsv.gz with q-values (collect-eqtl.py).
#
#   sbatch scripts/05-eqtl/collect-eqtl.sbatch <cohort> <ref> [method=star] [variants=matched]
set -euo pipefail
COHORT=${1:?cohort}; REF=${2:?ref}; METHOD=${3:-star}; VARIANTS=${4:-matched}
ROOT=/storage/home/yusuf.maulana/research/t2t-eqtl
export TMPDIR=/scratch/$USER/t2t-eqtl/tmp/$SLURM_JOB_ID; mkdir -p $TMPDIR; trap 'rm -rf $TMPDIR' EXIT
echo "start $(date -u +%FT%TZ) $COHORT $REF $METHOD $VARIANTS"
cd $ROOT
/usr/bin/time -v pixi run python scripts/05-eqtl/collect-eqtl.py --cohort $COHORT --ref $REF --method $METHOD --variants $VARIANTS
echo "end $(date -u +%FT%TZ)"
