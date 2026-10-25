#!/bin/bash
#SBATCH --job-name=t2t-tqtl
#SBATCH --partition=cpu
#SBATCH --array=0-22%6
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=/storage/projects/yusuf.maulana/t2t-eqtl/logs/05-eqtl/%x_%A_%a.out
#SBATCH --error=/storage/projects/yusuf.maulana/t2t-eqtl/logs/05-eqtl/%x_%A_%a.err
# Step 5: cis-eQTL mapping, one array task per chromosome (0-21 = chr1-22, 22 = chrX), at most 6 at a time.
#
#   sbatch scripts/05-eqtl/run-tensorqtl.sbatch <cohort> <ref> [method=star] [variants=matched] [nperm=10000]
#   then: pixi run python scripts/05-eqtl/collect-eqtl.py --cohort <cohort> --ref <ref> --method <m> --variants <v>
#
# Needs phenotypes.bed.gz and covariates.tsv (prepare-phenotypes.py, make-covariates.py) and the cohort's genotypes
# (prepare-genotypes.sbatch). CPU only (torch threads = CPUs); the node's GPU is used by the t2t-atlas runs.
set -euo pipefail
COHORT=${1:?cohort}; REF=${2:?ref}; METHOD=${3:-star}; VARIANTS=${4:-matched}; NPERM=${5:-10000}
ROOT=/storage/home/yusuf.maulana/research/t2t-eqtl
CHROMS=($(printf 'chr%s ' {1..22} X))
CHROM=${CHROMS[$SLURM_ARRAY_TASK_ID]}
export TMPDIR=/scratch/$USER/t2t-eqtl/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}; mkdir -p $TMPDIR; trap 'rm -rf $TMPDIR' EXIT
echo "start $(date -u +%FT%TZ) $COHORT $REF $METHOD $VARIANTS $CHROM"
/usr/bin/time -v pixi run --manifest-path $ROOT/pixi.toml python $ROOT/scripts/05-eqtl/run-tensorqtl.py \
  --cohort $COHORT --ref $REF --method $METHOD --variants $VARIANTS --chrom $CHROM --nperm $NPERM \
  --threads $SLURM_CPUS_PER_TASK
echo "end $(date -u +%FT%TZ)"
