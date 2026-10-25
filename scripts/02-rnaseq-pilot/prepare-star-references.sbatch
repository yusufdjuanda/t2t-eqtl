#!/bin/bash
#SBATCH --job-name=t2t-starref
#SBATCH --partition=cpu
#SBATCH --array=0-1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/storage/projects/yusuf.maulana/t2t-eqtl/logs/02-rnaseq-pilot/%x_%A_%a.out
#SBATCH --error=/storage/projects/yusuf.maulana/t2t-eqtl/logs/02-rnaseq-pilot/%x_%A_%a.err
# Step 2: build matched STAR references for RNA-seq, one per array task (0 = chm13v2.0, 1 = grch38).
#
#   sbatch scripts/02-rnaseq-pilot/prepare-star-references.sbatch
#
# Both use NCBI RefSeq annotation release RS_2025_08 (genes matched by NCBI GeneID). The T2T genome is
# changed only to remove two asymmetries with GRCh38_noalt (see docs/t2t-eqtl.md, step 2):
#   - chrY PARs hard-masked (N), as in GRCh38_noalt and the UCSC chm13v2.0_maskedY analysis set
#     (0-based: NC_060948.1 0-2458320 and 62122809-62460029)
#   - chrM and chrEBV appended, copied from GRCh38_noalt.fa; chrM genes copied from the GRCh38 GTF
# GRCh38 is used as is (links). Output in data/star-references/<id>/:
#   genome.fa(.fai)  annotation.gtf  star-2.7.11b/  transcripts.fa  genes.tsv  tx2gene.tsv  CHECKSUMS.md5
set -euo pipefail
ROOT=/storage/home/yusuf.maulana/research/t2t-eqtl
REFS=$ROOT/data/references
T2T_FA=$REFS/chm13v2.0-ncbi/genome/T2T-CHM13v2.0.fa
T2T_GTF=$REFS/chm13v2.0-ncbi/annotations/ncbi/T2T-CHM13v2.0.ncbiRefSeq.gtf
G38_FA=$REFS/grch38-noalt/genome/GRCh38_noalt.fa
G38_GTF=$REFS/grch38-noalt/annotations/ncbi/GRCh38_noalt.ncbiRefSeq.gtf
IDS=(chm13v2.0 grch38)
ID=${IDS[$SLURM_ARRAY_TASK_ID]}
OUT=$ROOT/data/star-references/$ID
module load samtools/1.23.1 bedtools/2.31.1 star/2.7.11b
export TMPDIR=/scratch/$USER/t2t-eqtl/tmp/${SLURM_JOB_ID}; mkdir -p $TMPDIR; trap 'rm -rf $TMPDIR' EXIT
mkdir -p $OUT
echo "start $(date -u +%FT%TZ) $ID"

if [[ $ID == chm13v2.0 ]]; then
  printf 'NC_060948.1\t0\t2458320\nNC_060948.1\t62122809\t62460029\n' > $TMPDIR/chrY-par.bed
  bedtools maskfasta -fi $T2T_FA -bed $TMPDIR/chrY-par.bed -fo $TMPDIR/masked.fa
  { cat $TMPDIR/masked.fa; samtools faidx $G38_FA chrM chrEBV; } > $OUT/genome.fa
  { grep -v '^#' $T2T_GTF; awk -F'\t' '$1=="chrM"' $G38_GTF; } > $OUT/annotation.gtf
else
  ln -sfn $G38_FA $OUT/genome.fa
  ln -sfn $G38_GTF $OUT/annotation.gtf
fi
samtools faidx $OUT/genome.fa

# gene_id  gene_ncbi  symbol  biotype  chrom  start  end  strand   (gene_ncbi = NCBI GeneID, the cross-reference key)
awk -F'\t' 'BEGIN{OFS="\t"; print "gene_id","gene_ncbi","symbol","biotype","chrom","start","end","strand"}
  $3=="gene" { id=$9; sub(/.*gene_id "/,"",id); sub(/".*/,"",id)
    g=""; if (match($9,/GeneID:[0-9]+/)) g=substr($9,RSTART+7,RLENGTH-7)
    s=id; if (match($9,/; gene "[^"]+"/)) s=substr($9,RSTART+8,RLENGTH-9)
    b=$9; sub(/.*gene_biotype "/,"",b); sub(/".*/,"",b)
    print id,g,s,b,$1,$4,$5,$7 }' $OUT/annotation.gtf > $OUT/genes.tsv
awk -F'\t' 'BEGIN{OFS="\t"; print "transcript_id","gene_id"}
  $3=="transcript" { t=$9; sub(/.*transcript_id "/,"",t); sub(/".*/,"",t)
    g=$9; sub(/.*gene_id "/,"",g); sub(/".*/,"",g); print t,g }' $OUT/annotation.gtf > $OUT/tx2gene.tsv
# gffread rejects RefSeq "gene" records (transcript_id ""), so give it transcript/exon/CDS records only
awk -F'\t' '$3!="gene"' $OUT/annotation.gtf > $TMPDIR/no-gene.gtf
pixi run --manifest-path $ROOT/pixi.toml gffread -w $OUT/transcripts.fa -g $OUT/genome.fa $TMPDIR/no-gene.gtf

# --sjdbOverhang 75: Geuvadis reads are 2 x 76 bp
rm -rf $OUT/star-2.7.11b
/usr/bin/time -v STAR --runMode genomeGenerate --runThreadN $SLURM_CPUS_PER_TASK \
  --genomeDir $OUT/star-2.7.11b --genomeFastaFiles $OUT/genome.fa \
  --sjdbGTFfile $OUT/annotation.gtf --sjdbOverhang 75 \
  --limitGenomeGenerateRAM 60000000000 --outTmpDir $TMPDIR/star --outFileNamePrefix $OUT/star-2.7.11b/

( cd $OUT && md5sum genome.fa genome.fa.fai annotation.gtf transcripts.fa genes.tsv tx2gene.tsv star-2.7.11b/* > CHECKSUMS.md5 )
echo "genes: $(($(wc -l < $OUT/genes.tsv)-1))  transcripts: $(grep -c '>' $OUT/transcripts.fa)  contigs: $(wc -l < $OUT/genome.fa.fai)"
echo "end $(date -u +%FT%TZ) $ID"
