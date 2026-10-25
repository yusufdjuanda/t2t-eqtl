#!/usr/bin/env bash
# Step 3: submit a whole cohort as a chain of Slurm jobs, batch by batch.
#
#   bash scripts/03-expression/submit-cohort.sh <cohort> [first_batch] [last_batch]
#   PREV_FETCH=<job> PREV_RUN=<job> bash … <cohort> <first> [last]   (attach to batch first-1 already queued/running)
#
# For each batch k: fetch_k (fetch-sra-batch.sbatch) then run_k (run-batch.sbatch).
#   fetch_k starts after fetch_{k-1} and run_{k-2} ended   -> at most ~2 batches of .sra on scratch
#   run_k   starts after fetch_k and run_{k-1} ended        -> one Nextflow head at a time
# All dependencies are afterany, so a failed download or batch never blocks the rest: Nextflow skips runs without a
# validated .sra (with a warning), and both scripts skip completed runs. Rerun a batch to pick up what is left:
#   sbatch scripts/03-expression/run-batch.sbatch <cohort> <k>   (resumes). Job IDs go to logs/03-expression/.
set -euo pipefail
COHORT=${1:?cohort}
ROOT=/storage/home/yusuf.maulana/research/t2t-eqtl
LAST_IN_SHEET=$(awk -F, 'NR==1 {for (i=1;i<=NF;i++) c[$i]=i; next} {print $c["batch"]}' $ROOT/config/samplesheets/$COHORT.csv | sort -n | tail -1)
FIRST=${2:-1}; LAST=${3:-$LAST_IN_SHEET}
LOG=$ROOT/logs/03-expression/submitted-$COHORT-$(date -u +%Y%m%dT%H%M%SZ).tsv
declare -A F R
[[ -n ${PREV_FETCH:-} ]] && F[$((FIRST-1))]=$PREV_FETCH
[[ -n ${PREV_RUN:-} ]] && R[$((FIRST-1))]=$PREV_RUN
printf 'batch\tfetch_job\trun_job\n' > $LOG
for ((k=FIRST; k<=LAST; k++)); do
  dep=""
  [[ -n ${F[$((k-1))]:-} ]] && dep="afterany:${F[$((k-1))]}"
  [[ -n ${R[$((k-2))]:-} ]] && dep="${dep:+$dep,}afterany:${R[$((k-2))]}"
  F[$k]=$(sbatch --parsable ${dep:+--dependency=$dep} $ROOT/scripts/03-expression/fetch-sra-batch.sbatch $COHORT $k)
  dep="afterany:${F[$k]}"
  [[ -n ${R[$((k-1))]:-} ]] && dep="$dep,afterany:${R[$((k-1))]}"
  R[$k]=$(sbatch --parsable --dependency=$dep $ROOT/scripts/03-expression/run-batch.sbatch $COHORT $k)
  printf '%s\t%s\t%s\n' $k ${F[$k]} ${R[$k]} | tee -a $LOG
done
echo "submitted batches $FIRST-$LAST of $COHORT; job list in $LOG"
