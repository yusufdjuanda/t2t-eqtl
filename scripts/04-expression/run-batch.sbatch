#!/bin/bash
#SBATCH --job-name=t2t-nf3
# t2t-partition: long
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/04-expression/%x_%j.out
#SBATCH --error=logs/04-expression/%x_%j.err
# Stage 04 (expression), part 2: Nextflow head job for one batch; then delete the .sra of every run whose results are complete.
#
#   scripts/submit.sh scripts/04-expression/run-batch.sbatch <cohort> <batch> [samplesheet]
#   scripts/submit.sh scripts/04-expression/run-batch.sbatch mage 1 test        (2-run test: profile test)
#
# Each batch has its own launch folder (work/04-expression/launch/<cohort>-b<batch>/), so rerunning the same
# command resumes that batch: finished tasks are reused from the shared work folder, and runs with
# complete.txt are skipped altogether. Trace: logs/04-expression/trace-<cohort>-b<batch>-<jobid>.tsv.
set -euo pipefail
COHORT=${1:?cohort}; BATCH=${2:?batch}; SHEET=${3:-}
ROOT=${T2T_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && git rev-parse --show-toplevel)}}
source "$ROOT/scripts/lib/common.sh"
PROFILE=slurm
if [[ $SHEET == test ]]; then PROFILE=slurm,test; SHEET=$ROOT/config/samplesheets/test-mage-2.csv; fi
SHEET=${SHEET:-$ROOT/config/samplesheets/$COHORT.csv}
LAUNCH=$ROOT/work/04-expression/launch/$COHORT-b$BATCH
SRA=$ROOT/work/04-expression/sra/$COHORT
RES=$ROOT/results/04-expression/$COHORT
t2t_env bio   # STAR, salmon, samtools, ... from pixi.toml (published runs: cluster modules of the same versions)
t2t_tmpdir "$SLURM_JOB_ID"
export NXF_OPTS='-Xms1g -Xmx4g'
mkdir -p $LAUNCH; cd $LAUNCH
echo "start $(date -u +%FT%TZ) $COHORT batch $BATCH"

nextflow run $ROOT/scripts/04-expression/nextflow/main.nf -profile $PROFILE \
  --input $SHEET --batch $BATCH \
  -work-dir $ROOT/work/04-expression/nf-work/$COHORT -resume \
  -with-trace $ROOT/logs/04-expression/trace-$COHORT-b$BATCH-$SLURM_JOB_ID.tsv

# delete .sra only for runs whose results are complete
n=0
for d in $SRA/*/; do
  run=$(basename $d)
  if [[ -f $RES/$run/complete.txt ]]; then rm -rf $d; n=$((n+1)); fi
done
echo "deleted .sra of $n completed runs; $(du -sh $SRA 2>/dev/null | cut -f1) left in $SRA"
echo "end $(date -u +%FT%TZ)"
