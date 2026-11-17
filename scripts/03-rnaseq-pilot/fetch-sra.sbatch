#!/bin/bash
#SBATCH --job-name=t2t-sra
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=logs/03-rnaseq-pilot/%x_%j.out
#SBATCH --error=logs/03-rnaseq-pilot/%x_%j.err
# Step 2: get Geuvadis reads from NCBI SRA and write them as gzipped paired FASTQs.
#
#   scripts/submit.sh scripts/03-rnaseq-pilot/fetch-sra.sbatch [sample sheet]
#
# Replaces fetch-fastq.sbatch (ENA), which gave one connection at ~0.75 MB/s on 2026-10-17; NCBI gave ~12 MB/s.
# NCBI serves .sra files from AWS S3, where connections from Karakoram are often reset, so prefetch is retried
# (it resumes partial downloads). Per run, 3 runs in parallel (4 threads each):
#   prefetch (resume) -> vdb-validate (NCBI checksums) -> fasterq-dump --split-files -> pigz
#   -> read pairs must equal read_pairs in the sheet (ENA's count) -> FASTQs read-only -> .sra deleted
# Output: work/03-rnaseq-pilot/fastq/<run>/<run>_{1,2}.fastq.gz. A run with read-only FASTQs is skipped on rerun.
set -euo pipefail
ROOT=${T2T_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && git rev-parse --show-toplevel)}}
source "$ROOT/scripts/lib/common.sh"
SHEET=${1:-$ROOT/config/samples/pilot-5.tsv}
OUT=$ROOT/work/03-rnaseq-pilot/fastq
SRA=$ROOT/work/03-rnaseq-pilot/sra
t2t_env bio   # STAR, salmon, samtools, ... from pixi.toml (published runs: cluster modules of the same versions)
t2t_tmpdir "$SLURM_JOB_ID"
mkdir -p $OUT $SRA
echo "start $(date -u +%FT%TZ)"

fetch_one() {  # run expected_pairs
  local run=$1 pairs=$2 d=$OUT/$1 n
  if [[ -f $d/${run}_2.fastq.gz && ! -w $d/${run}_2.fastq.gz ]]; then echo "$run already done"; return 0; fi
  mkdir -p $d
  for i in $(seq 1 100); do
    prefetch --max-size 50G -O $SRA $run >> $TMPDIR/$run.prefetch.log 2>&1 && break
    echo "$(date -u +%T) $run prefetch attempt $i interrupted at $(du -sb $SRA/$run 2>/dev/null | cut -f1) bytes" >&2; sleep 10
  done
  vdb-validate $SRA/$run/$run.sra 2>&1 | tail -n 3 >&2
  fasterq-dump --split-files --threads 4 --temp $TMPDIR -O $d $SRA/$run/$run.sra
  pigz -p 4 $d/${run}_1.fastq $d/${run}_2.fastq
  n=$(( $(pigz -dc -p 4 $d/${run}_1.fastq.gz | wc -l) / 4 ))
  if [[ $n -ne $pairs ]]; then echo "$run read pairs $n != expected $pairs" >&2; return 1; fi
  chmod a-w $d/${run}_1.fastq.gz $d/${run}_2.fastq.gz
  rm -rf $SRA/$run
  echo "$(date -u +%T) $run ok: $n read pairs (matches ENA), $(du -cb $d/*.fastq.gz | tail -1 | cut -f1) bytes"
}
export -f fetch_one; export OUT SRA TMPDIR

awk -F'\t' 'NR>1 {print $1, $6}' $SHEET \
  | /usr/bin/time -v xargs -P 3 -L1 bash -c 'set -euo pipefail; fetch_one "$@"' _
du -sb $OUT
echo "end $(date -u +%FT%TZ)"
