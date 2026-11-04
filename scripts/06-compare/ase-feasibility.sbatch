#!/bin/bash
#SBATCH --job-name=t2t-asefeas
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=/storage/projects/yusuf.maulana/t2t-eqtl/logs/06-compare/%x_%j.out
#SBATCH --error=/storage/projects/yusuf.maulana/t2t-eqtl/logs/06-compare/%x_%j.err
# Step 6: feasibility of allele-specific expression at the changed genes (ase-feasibility.py).
#
#   sbatch scripts/06-compare/ase-feasibility.sbatch <cohort>
set -euo pipefail
COHORT=${1:?cohort}
ROOT=/storage/home/yusuf.maulana/research/t2t-eqtl
export TMPDIR=/scratch/$USER/t2t-eqtl/tmp/$SLURM_JOB_ID; mkdir -p $TMPDIR; trap 'rm -rf $TMPDIR' EXIT
module load minimap2/2.31
echo "start $(date -u +%FT%TZ) $COHORT"
cd $ROOT
/usr/bin/time -v pixi run python scripts/06-compare/ase-feasibility.py --cohort $COHORT --threads $SLURM_CPUS_PER_TASK
echo "end $(date -u +%FT%TZ)"
