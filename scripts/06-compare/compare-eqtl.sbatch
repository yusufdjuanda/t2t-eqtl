#!/bin/bash
#SBATCH --job-name=t2t-cmpeqtl
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=logs/06-compare/%x_%j.out
#SBATCH --error=logs/06-compare/%x_%j.err
# Step 6: compare the eQTL results of the three references for one cohort, quantification method and variant set.
#
#   scripts/submit.sh scripts/06-compare/compare-eqtl.sbatch <cohort> [method=star] [variants=matched] [skip-genes]
#
# Runs compare-eqtl.py (eGene presence, lead-variant replication -> genes.tsv, summary.tsv) and then
# compare-effects.py (same-variant effect sizes -> effects-*.tsv). With a 4th argument "skip-genes" only the second
# runs (genes.tsv must exist). Needs collect-eqtl.py to have run for all three references.
set -euo pipefail
COHORT=${1:?cohort}; METHOD=${2:-star}; VARIANTS=${3:-matched}; SKIP=${4:-}
ROOT=${T2T_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && git rev-parse --show-toplevel)}}
source "$ROOT/scripts/lib/common.sh"
t2t_tmpdir "$SLURM_JOB_ID"
echo "start $(date -u +%FT%TZ) $COHORT $METHOD $VARIANTS"
cd $ROOT
if [[ $SKIP != skip-genes ]]; then
  /usr/bin/time -v pixi run python scripts/06-compare/compare-eqtl.py --cohort $COHORT --method $METHOD --variants $VARIANTS
fi
/usr/bin/time -v pixi run python scripts/06-compare/compare-effects.py --cohort $COHORT --method $METHOD --variants $VARIANTS
echo "end $(date -u +%FT%TZ)"
