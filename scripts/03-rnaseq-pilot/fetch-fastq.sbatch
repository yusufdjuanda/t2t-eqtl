#!/bin/bash
#SBATCH --job-name=t2t-fetch
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=/storage/projects/yusuf.maulana/t2t-eqtl/logs/03-rnaseq-pilot/%x_%j.out
#SBATCH --error=/storage/projects/yusuf.maulana/t2t-eqtl/logs/03-rnaseq-pilot/%x_%j.err
# Step 2: download Geuvadis FASTQs from ENA and check their md5 against ENA.
#
#   sbatch scripts/03-rnaseq-pilot/fetch-fastq.sbatch [sample sheet]
#
# Sample sheet (default config/samples/pilot-5.tsv): run sample ancestry sex lab read_pairs url_1 url_2 md5_1 md5_2.
# Files go to work/03-rnaseq-pilot/fastq/<run>/ (scratch), one stream per CPU (default 2), resumed on failure.
# ENA refused or reset most connections when 10 streams were opened at once (2026-10-17), so keep this low.
# A verified file is made read-only; rerunning skips it. Exits non-zero if any md5 fails.
set -euo pipefail
ROOT=/storage/home/yusuf.maulana/research/t2t-eqtl
SHEET=${1:-$ROOT/config/samples/pilot-5.tsv}
OUT=$ROOT/work/03-rnaseq-pilot/fastq
export TMPDIR=/scratch/$USER/t2t-eqtl/tmp/$SLURM_JOB_ID; mkdir -p $TMPDIR; trap 'rm -rf $TMPDIR' EXIT
echo "start $(date -u +%FT%TZ)"

fetch_one() {  # run url md5
  local run=$1 url=$2 md5=$3 f
  f=$OUT/$run/$(basename $url)
  mkdir -p $OUT/$run
  if [[ -f $f && ! -w $f ]]; then echo "$f already verified"; return 0; fi
  for i in $(seq 1 60); do
    curl -sSf -C - --retry 3 --retry-all-errors --retry-delay 30 -o $f $url && break
    echo "$(date -u +%T) $f attempt $i failed at $(stat -c %s $f 2>/dev/null || echo 0) bytes" >&2; sleep 60
  done
  if [[ $(md5sum < $f | cut -d' ' -f1) == $md5 ]]; then chmod a-w $f; echo "$f md5 ok"
  else echo "$f md5 MISMATCH" >&2; return 1; fi
}
export -f fetch_one; export OUT

awk -F'\t' 'NR>1 {print $1, $7, $9; print $1, $8, $10}' $SHEET \
  | /usr/bin/time -v xargs -P $SLURM_CPUS_PER_TASK -L1 bash -c 'fetch_one "$@"' _
du -sb $OUT
echo "end $(date -u +%FT%TZ)"
