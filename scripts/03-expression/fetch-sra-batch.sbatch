#!/bin/bash
#SBATCH --job-name=t2t-fetch3
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=6
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/storage/projects/yusuf.maulana/t2t-eqtl/logs/03-expression/%x_%j.out
#SBATCH --error=/storage/projects/yusuf.maulana/t2t-eqtl/logs/03-expression/%x_%j.err
# Step 3, part 1: download one batch of runs from NCBI SRA (the .sra files only; Nextflow converts them to FASTQ).
#
#   sbatch scripts/03-expression/fetch-sra-batch.sbatch <cohort> <batch> [samplesheet]      cohort = mage | geuvadis
#
# Runs = rows of config/samplesheets/<cohort>.csv with that batch number (60 runs, ~160 GB of .sra); the .sra goes
# to the samplesheet's `sra` path (work/03-expression/sra/<cohort>/<run>/<run>.sra).
# Per run, 3 in parallel: prefetch (retried; resumes after the AWS S3 connection resets seen from Karakoram)
# -> vdb-validate (NCBI checksums) -> marker <run>.validated. Skipped: runs already validated, and runs whose
# results are complete (results/03-expression/<cohort>/<run>/complete.txt). Before each download it waits while
# scratch has less than MIN_FREE_GB (default 400) free. Exits non-zero if any run fails.
# Output: work/03-expression/sra/<cohort>/<run>/<run>.sra (scratch). run-batch.sbatch deletes each .sra once that
# run's results are complete, so at most about two batches are on disk at a time (submit-cohort.sh).
set -euo pipefail
COHORT=${1:?cohort}; BATCH=${2:?batch}
ROOT=/storage/home/yusuf.maulana/research/t2t-eqtl
SHEET=${3:-$ROOT/config/samplesheets/$COHORT.csv}
SRA=$ROOT/work/03-expression/sra/$COHORT
MIN_FREE_GB=${MIN_FREE_GB:-400}   # do not download while scratch has less free space than this
RES=$ROOT/results/03-expression/$COHORT
module load sra-tools/3.4.1
export TMPDIR=/scratch/$USER/t2t-eqtl/tmp/$SLURM_JOB_ID; mkdir -p $TMPDIR; trap 'rm -rf $TMPDIR' EXIT
mkdir -p $SRA
echo "start $(date -u +%FT%TZ) $COHORT batch $BATCH"

fetch_one() {  # run
  local run=$1
  if [[ -f $RES/$run/complete.txt ]]; then echo "$run results complete, skipped"; return 0; fi
  if [[ -f $SRA/$run/$run.validated ]]; then echo "$run already validated"; return 0; fi
  # scratch is shared (it reached 99% on 2026-10-23): wait while less than MIN_FREE_GB is free before downloading
  while (( $(df --output=avail -BG $SRA | tail -1 | tr -dc 0-9) < MIN_FREE_GB )); do
    echo "$(date -u +%T) $run waiting: scratch has $(df --output=avail -h $SRA | tail -1 | tr -d ' ') free (< ${MIN_FREE_GB}G)" >&2
    sleep 600
  done
  # prefetch can loop for hours without failing (2026-10-17: stuck at ~120 MB, restarting its own partial file),
  # so each attempt gets a time limit; after a timed-out attempt the partial files are removed and it starts fresh.
  for i in $(seq 1 100); do
    rc=0; timeout 30m prefetch --max-size 100G -O $SRA $run >> $TMPDIR/$run.prefetch.log 2>&1 || rc=$?
    [[ $rc -eq 0 ]] && break
    echo "$(date -u +%T) $run prefetch attempt $i exit $rc at $(du -sb $SRA/$run 2>/dev/null | cut -f1) bytes" >&2
    if [[ $rc -eq 124 ]]; then
      rm -f -- "${SRA:?}/${run:?}/${run:?}.sra.tmp" "${SRA:?}/${run:?}/${run:?}.sra.prf" "${SRA:?}/${run:?}/${run:?}.sra.lock"
    fi
    sleep 10
  done
  vdb-validate $SRA/$run/$run.sra > $TMPDIR/$run.validate.log 2>&1 || { cat $TMPDIR/$run.validate.log >&2; return 1; }
  date -u +%FT%TZ > $SRA/$run/$run.validated
  echo "$(date -u +%T) $run validated, $(stat -c %s $SRA/$run/$run.sra) bytes"
}
export -f fetch_one; export SRA RES TMPDIR MIN_FREE_GB

# column positions by name, so the samplesheet column order may change
awk -F, -v b=$BATCH 'NR==1 {for (i=1;i<=NF;i++) c[$i]=i; next} $c["batch"]==b {print $c["run"]}' $SHEET \
  | /usr/bin/time -v xargs -P 3 -I{} bash -c 'set -euo pipefail; fetch_one {}'
echo "batch $BATCH: $(ls $SRA/*/*.validated 2>/dev/null | wc -l) validated .sra on disk, $(du -sh $SRA | cut -f1)"
echo "end $(date -u +%FT%TZ)"
