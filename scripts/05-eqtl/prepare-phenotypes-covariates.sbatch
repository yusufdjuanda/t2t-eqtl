#!/bin/bash
#SBATCH --job-name=t2t-phenocov
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=logs/05-eqtl/%x_%j.out
#SBATCH --error=logs/05-eqtl/%x_%j.err
# Step 5: phenotype BED and covariates for one cohort and quantification method, on all three references.
#
#   scripts/submit.sh scripts/05-eqtl/prepare-phenotypes-covariates.sbatch <cohort> [method=star]
#
# Runs prepare-phenotypes.py then make-covariates.py per reference (genotype PCs are reused from the cohort cache).
set -euo pipefail
COHORT=${1:?cohort}; METHOD=${2:-star}
ROOT=${T2T_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && git rev-parse --show-toplevel)}}
source "$ROOT/scripts/lib/common.sh"
t2t_tmpdir "$SLURM_JOB_ID"
echo "start $(date -u +%FT%TZ) $COHORT $METHOD"
cd $ROOT
for r in chm13v2.0 grch38 grch38-ffx; do
  /usr/bin/time -v pixi run python scripts/05-eqtl/prepare-phenotypes.py --cohort $COHORT --ref $r --method $METHOD
  pixi run python scripts/05-eqtl/make-covariates.py --cohort $COHORT --ref $r --method $METHOD
done
echo "end $(date -u +%FT%TZ)"
