#!/bin/bash
#SBATCH --job-name=t2t-rnaaln
#SBATCH --partition=cpu
#SBATCH --array=0-9%2
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=/storage/projects/yusuf.maulana/t2t-eqtl/logs/02-rnaseq-pilot/%x_%A_%a.out
#SBATCH --error=/storage/projects/yusuf.maulana/t2t-eqtl/logs/02-rnaseq-pilot/%x_%A_%a.err
# Step 2: align one Geuvadis run to one reference with STAR and quantify it two ways.
# Array task i = sample row i/2 of the sheet x reference i%2 (0 = chm13v2.0, 1 = grch38); 5 samples -> 0-9.
#
#   sbatch scripts/02-rnaseq-pilot/align-quantify.sbatch [sample sheet]
#
# Inputs: verified FASTQs from fetch-fastq.sbatch; data/star-references/<ref>/ from prepare-star-references.sbatch.
# Same STAR version, index parameters and alignment parameters (ENCODE long-RNA settings) on both references.
#   STAR GeneCounts  uniquely mapped read pairs per gene (multi-mappers dropped). Geuvadis is unstranded: column 2.
#   salmon (alignment mode, on STAR's transcriptome BAM)  multi-mapping reads shared between transcripts by EM.
# Output: results/02-rnaseq-pilot/<ref>/<run>/  Log.final.out  ReadsPerGene.out.tab  SJ.out.tab  idxstats.tsv  salmon/
#         work/02-rnaseq-pilot/bam/<ref>/<run>.bam(.bai)  genome BAM on scratch, for inspecting loci
set -euo pipefail
ROOT=/storage/home/yusuf.maulana/research/t2t-eqtl
SHEET=${1:-$ROOT/config/samples/pilot-5.tsv}
REFS_LIST=(chm13v2.0 grch38)
REF=${REFS_LIST[$((SLURM_ARRAY_TASK_ID % 2))]}
RUN=$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID / 2 + 2)) 'NR==n {print $1}' $SHEET)
[[ -n $RUN ]] || { echo "no sample for task $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
SREF=$ROOT/data/star-references/$REF
FQ=$ROOT/work/02-rnaseq-pilot/fastq/$RUN
OUT=$ROOT/results/02-rnaseq-pilot/$REF/$RUN
BAMDIR=$ROOT/work/02-rnaseq-pilot/bam/$REF
module load star/2.7.11b salmon/1.11.4 samtools/1.23.1
export TMPDIR=/scratch/$USER/t2t-eqtl/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}; mkdir -p $TMPDIR; trap 'rm -rf $TMPDIR' EXIT
mkdir -p $OUT $BAMDIR
echo "start $(date -u +%FT%TZ) $RUN $REF"

/usr/bin/time -v STAR --runThreadN $SLURM_CPUS_PER_TASK --genomeDir $SREF/star-2.7.11b \
  --readFilesIn $FQ/${RUN}_1.fastq.gz $FQ/${RUN}_2.fastq.gz --readFilesCommand zcat \
  --outFileNamePrefix $TMPDIR/ --outTmpDir $TMPDIR/star \
  --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 \
  --outFilterMismatchNoverReadLmax 0.04 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
  --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
  --outSAMattributes NH HI AS NM MD --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 20000000000 \
  --quantMode TranscriptomeSAM GeneCounts

/usr/bin/time -v salmon quant -p $SLURM_CPUS_PER_TASK -l A -t $SREF/transcripts.fa \
  -a $TMPDIR/Aligned.toTranscriptome.out.bam -o $OUT/salmon

mv $TMPDIR/Aligned.sortedByCoord.out.bam $BAMDIR/$RUN.bam
samtools index -@ $SLURM_CPUS_PER_TASK $BAMDIR/$RUN.bam
samtools idxstats $BAMDIR/$RUN.bam > $OUT/idxstats.tsv
cp $TMPDIR/Log.final.out $TMPDIR/Log.out $TMPDIR/ReadsPerGene.out.tab $TMPDIR/SJ.out.tab $OUT/
du -sb $FQ $BAMDIR/$RUN.bam $OUT
echo "end $(date -u +%FT%TZ) $RUN $REF"
