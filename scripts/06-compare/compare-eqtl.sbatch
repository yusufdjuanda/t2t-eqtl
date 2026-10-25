#!/bin/bash
#SBATCH --job-name=t2t-cmpeqtl
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/storage/projects/yusuf.maulana/t2t-eqtl/logs/06-compare/%x_%j.out
#SBATCH --error=/storage/projects/yusuf.maulana/t2t-eqtl/logs/06-compare/%x_%j.err
# Step 6: compare the eQTL results of the three references for one cohort, quantification method and variant set.
#
#   sbatch scripts/06-compare/compare-eqtl.sbatch <cohort> [method=star] [variants=matched] [skip-genes]
#
# Runs compare-eqtl.py (eGene presence, lead-variant replication -> genes.tsv, summary.tsv) and then
# compare-effects.py (same-variant effect sizes -> effects-*.tsv). With a 4th argument "skip-genes" only the second
# runs (genes.tsv must exist). Needs collect-eqtl.py to have run for all three references.
set -euo pipefail
COHORT=${1:?cohort}; METHOD=${2:-star}; VARIANTS=${3:-matched}; SKIP=${4:-}
ROOT=/storage/home/yusuf.maulana/research/t2t-eqtl
export TMPDIR=/scratch/$USER/t2t-eqtl/tmp/$SLURM_JOB_ID; mkdir -p $TMPDIR; trap 'rm -rf $TMPDIR' EXIT
echo "start $(date -u +%FT%TZ) $COHORT $METHOD $VARIANTS"
cd $ROOT
if [[ $SKIP != skip-genes ]]; then
  /usr/bin/time -v pixi run python scripts/06-compare/compare-eqtl.py --cohort $COHORT --method $METHOD --variants $VARIANTS
fi
/usr/bin/time -v pixi run python scripts/06-compare/compare-effects.py --cohort $COHORT --method $METHOD --variants $VARIANTS
echo "end $(date -u +%FT%TZ)"
