#!/bin/bash
#SBATCH --job-name=t2t-aserna
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=/storage/projects/yusuf.maulana/t2t-eqtl/logs/06-compare/%x_%A_%a.out
#SBATCH --error=/storage/projects/yusuf.maulana/t2t-eqtl/logs/06-compare/%x_%A_%a.err
# Step 6: RNA-seq realignment for allele-specific expression (ASE), one run per array task.
#
#   sbatch --array=1-<n>%10 scripts/06-compare/ase-realign.sbatch [runs.tsv]
#
# Runs: results/06-compare/ase/targets/runs.tsv (ase-targets.py), task i = data row i.
# Per run: prefetch from SRA (retried, as fetch-sra-batch.sbatch) -> vdb-validate -> fasterq-dump -> read-pair check ->
# for chm13v2.0 and grch38-ffx: STAR with exactly the stage-04 settings (QUANTIFY_RUN), no salmon -> keep the reads in
# targets/rna-<ref>.bed (all alignments, including multi-mappers) -> delete the genome BAM, FASTQ and .sra.
# Output: results/06-compare/ase/rna/<run>/<ref>/{targets.bam(.bai), Log.final.out}, complete.txt (runs that have it
# are skipped). Waits while scratch has less than MIN_FREE_GB (default 200) free.
set -euo pipefail
ROOT=/storage/home/yusuf.maulana/research/t2t-eqtl
RUNS=${1:-$ROOT/results/06-compare/ase/targets/runs.tsv}
MIN_FREE_GB=${MIN_FREE_GB:-200}
read -r RUN SAMPLE PAIRS RL _ < <(awk -F'\t' -v i=$SLURM_ARRAY_TASK_ID 'NR == i + 1' $RUNS)
OUT=$ROOT/results/06-compare/ase/rna/$RUN
[[ -f $OUT/complete.txt ]] && { echo "$RUN complete, skipped"; exit 0; }
module load sra-tools/3.4.1 star/2.7.11b samtools/1.23.1
export TMPDIR=/scratch/$USER/t2t-eqtl/tmp/ase-$SLURM_ARRAY_JOB_ID-$RUN; mkdir -p $TMPDIR; trap 'rm -rf $TMPDIR' EXIT
mkdir -p $OUT
echo "start $(date -u +%FT%TZ) $RUN $SAMPLE"

while (( $(df --output=avail -BG /scratch | tail -1 | tr -dc 0-9) < MIN_FREE_GB )); do
  echo "$(date -u +%T) waiting: scratch below ${MIN_FREE_GB}G free" >&2; sleep 600
done
for i in $(seq 1 100); do
  rc=0; timeout 30m prefetch --max-size 100G -O $TMPDIR $RUN >> $TMPDIR/prefetch.log 2>&1 || rc=$?
  [[ $rc -eq 0 ]] && break
  echo "$(date -u +%T) prefetch attempt $i exit $rc" >&2
  [[ $rc -eq 124 ]] && rm -f -- "${TMPDIR:?}/${RUN:?}/${RUN:?}.sra.tmp" "${TMPDIR:?}/${RUN:?}/${RUN:?}.sra.prf" "${TMPDIR:?}/${RUN:?}/${RUN:?}.sra.lock"
  sleep 10
done
vdb-validate $TMPDIR/$RUN/$RUN.sra > $TMPDIR/validate.log 2>&1 || { cat $TMPDIR/validate.log >&2; exit 1; }
echo "$(date -u +%T) downloaded $(stat -c %s $TMPDIR/$RUN/$RUN.sra) bytes"
fasterq-dump --split-files --threads $SLURM_CPUS_PER_TASK --temp $TMPDIR -O $TMPDIR/fq $TMPDIR/$RUN/$RUN.sra
rm -f -- "${TMPDIR:?}/${RUN:?}/${RUN:?}.sra"
n=$(( $(wc -l < $TMPDIR/fq/${RUN}_1.fastq) / 4 ))
[[ $n -eq $PAIRS ]] || { echo "$RUN: $n read pairs, expected $PAIRS" >&2; exit 1; }

IDX=star-2.7.11b-oh$((RL - 1))
for ref in chm13v2.0 grch38-ffx; do
  R=$ROOT/data/star-references/$ref
  O=$OUT/$ref; mkdir -p $O
  /usr/bin/time -v -o $O/time-star.txt STAR --runThreadN $SLURM_CPUS_PER_TASK --genomeDir $R/$IDX \
    --readFilesIn $TMPDIR/fq/${RUN}_1.fastq $TMPDIR/fq/${RUN}_2.fastq \
    --outFileNamePrefix $TMPDIR/$ref. --outTmpDir $TMPDIR/star-$ref \
    --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
    --outSAMattributes NH HI AS NM MD --outSAMunmapped Within \
    --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 20000000000
  samtools index -@ $SLURM_CPUS_PER_TASK $TMPDIR/$ref.Aligned.sortedByCoord.out.bam
  samtools view -@ $SLURM_CPUS_PER_TASK -b -M -L $ROOT/results/06-compare/ase/targets/rna-$ref.bed \
    -o $O/targets.bam $TMPDIR/$ref.Aligned.sortedByCoord.out.bam
  samtools index $O/targets.bam
  cp $TMPDIR/$ref.Log.final.out $O/Log.final.out
  rm -f -- "${TMPDIR:?}/${ref:?}".Aligned.*
done
echo "end $(date -u +%FT%TZ)"
date -u +%FT%TZ > $OUT/complete.txt
