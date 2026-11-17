#!/bin/bash
#SBATCH --job-name=t2t-tqtl
#SBATCH --array=0-22%6
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=logs/05-eqtl/%x_%A_%a.out
#SBATCH --error=logs/05-eqtl/%x_%A_%a.err
# Step 5: cis-eQTL mapping, one array task per chromosome (0-21 = chr1-22, 22 = chrX), at most 6 at a time.
#
#   scripts/submit.sh scripts/05-eqtl/run-tensorqtl.sbatch <cohort> <ref> [method=star] [variants=matched] [nperm=10000]
#   then: pixi run python scripts/05-eqtl/collect-eqtl.py --cohort <cohort> --ref <ref> --method <m> --variants <v>
#
# Needs phenotypes.bed.gz and covariates.tsv (prepare-phenotypes.py, make-covariates.py) and the cohort's genotypes
# (prepare-genotypes.sbatch). CPU only (torch threads = CPUs); the node's GPU is used by the t2t-atlas runs.
set -euo pipefail
COHORT=${1:?cohort}; REF=${2:?ref}; METHOD=${3:-star}; VARIANTS=${4:-matched}; NPERM=${5:-10000}
ROOT=${T2T_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && git rev-parse --show-toplevel)}}
source "$ROOT/scripts/lib/common.sh"
CHROMS=($(printf 'chr%s ' {1..22} X))
CHROM=${CHROMS[$SLURM_ARRAY_TASK_ID]}
t2t_tmpdir "${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "start $(date -u +%FT%TZ) $COHORT $REF $METHOD $VARIANTS $CHROM"
/usr/bin/time -v pixi run --manifest-path $ROOT/pixi.toml python $ROOT/scripts/05-eqtl/run-tensorqtl.py \
  --cohort $COHORT --ref $REF --method $METHOD --variants $VARIANTS --chrom $CHROM --nperm $NPERM \
  --threads $SLURM_CPUS_PER_TASK
echo "end $(date -u +%FT%TZ)"
