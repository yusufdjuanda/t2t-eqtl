#!/bin/bash
#SBATCH --job-name=t2t-collect
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/05-eqtl/%x_%j.out
#SBATCH --error=logs/05-eqtl/%x_%j.err
# Step 5: per-chromosome tensorQTL results -> cis-eqtl.tsv.gz with q-values (collect-eqtl.py).
#
#   scripts/submit.sh scripts/05-eqtl/collect-eqtl.sbatch <cohort> <ref> [method=star] [variants=matched]
set -euo pipefail
COHORT=${1:?cohort}; REF=${2:?ref}; METHOD=${3:-star}; VARIANTS=${4:-matched}
ROOT=${T2T_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && git rev-parse --show-toplevel)}}
source "$ROOT/scripts/lib/common.sh"
t2t_tmpdir "$SLURM_JOB_ID"
echo "start $(date -u +%FT%TZ) $COHORT $REF $METHOD $VARIANTS"
cd $ROOT
/usr/bin/time -v pixi run python scripts/05-eqtl/collect-eqtl.py --cohort $COHORT --ref $REF --method $METHOD --variants $VARIANTS
echo "end $(date -u +%FT%TZ)"
