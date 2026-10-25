#!/bin/bash
#SBATCH --job-name=t2t-geno
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=23
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=/storage/projects/yusuf.maulana/t2t-eqtl/logs/01-genotypes/%x_%j.out
#SBATCH --error=/storage/projects/yusuf.maulana/t2t-eqtl/logs/01-genotypes/%x_%j.err
# Step 1: subset the phased 1KGP T2T-CHM13v2.0 panel (SHAPEIT5 v1.1, 3,202 samples) to the
# Geuvadis samples, and to NA12878 for allele-specific binding, then count variants by class.
#
#   sbatch scripts/01-genotypes/subset-genotypes.sbatch [output folder] [sample list]
#   (step 3: ... results/01b-expression-individuals-genotypes-chm13 config/samples/expression-individuals-in-1kgp-chm13.txt)
#
# Per chromosome, in parallel (one bcftools pipeline per contig):
#   geuvadis/<chr>.bcf  the 451 Geuvadis samples in config/samples/geuvadis-in-1kgp-chm13.txt; multiallelic
#                       sites split (norm -m-), AC/AN/MAF recomputed in this subset, monomorphic sites dropped
#   na12878/<chr>.bcf   NA12878, sites where it carries a non-reference allele
#   counts/<chr>.tsv    variant counts by SYNTENIC flag x type x Geuvadis MAF bin, and NA12878 het counts
# SYNTENIC (INFO flag set by the panel authors) = in sequence syntenic with GRCh38; absent = T2T-only sequence
# (source: chm13v2-unique_to_hg38.bed). Output: summary.tsv, all chromosomes combined.
set -euo pipefail
ROOT=/storage/home/yusuf.maulana/research/t2t-eqtl
OUT=${1:-$ROOT/results/01-geuvadis-genotypes-chm13}
PANEL=$ROOT/data/external/1kgp-chm13v2.0-phased-shapeit5-v1.1/1KGP.CHM13v2.0.whole_genome.recalibrated.snp_indel.pass.phased.native_maps.3202.bcf.gz
SAMPLES=${2:-$ROOT/config/samples/geuvadis-in-1kgp-chm13.txt}
module load bcftools/1.23.1
export TMPDIR=/scratch/$USER/t2t-eqtl/tmp/$SLURM_JOB_ID; mkdir -p $TMPDIR; trap 'rm -rf $TMPDIR' EXIT
mkdir -p $OUT/geuvadis $OUT/na12878 $OUT/counts
echo "start $(date -u +%FT%TZ)"

one_chrom() {
  local c=$1
  bcftools view -S $SAMPLES -Ou -r $c $PANEL \
    | bcftools norm -m- -Ou \
    | bcftools +fill-tags -Ou -- -t AC,AN,MAF \
    | bcftools view -c1:minor -Ob -W=csi -o $OUT/geuvadis/$c.bcf
  bcftools view -s NA12878 -r $c -c1 -Ob -W=csi -o $OUT/na12878/$c.bcf $PANEL
  # class  type  maf_bin  n   (class: syntenic | t2t_only)
  bcftools query -f '%INFO/SYNTENIC\t%TYPE\t%INFO/MAF\n' $OUT/geuvadis/$c.bcf \
    | awk -v c=$c 'BEGIN{OFS="\t"} {cls=($1=="1")?"syntenic":"t2t_only"; t=($2=="SNP")?"snp":"indel";
        b=($3>=0.05)?"maf>=0.05":($3>=0.01)?"maf0.01-0.05":"maf<0.01"; n[cls"\t"t"\t"b]++}
        END{for(k in n) print c,"geuvadis",k,n[k]}' > $OUT/counts/$c.tsv
  bcftools query -f '%INFO/SYNTENIC\t[%GT]\n' $OUT/na12878/$c.bcf \
    | awk -v c=$c 'BEGIN{OFS="\t"} {cls=($1=="1")?"syntenic":"t2t_only"; split($2,g,/[|\/]/);
        z=(g[1]!=g[2])?"het":"hom_alt"; n[cls"\t"z]++}
        END{for(k in n) print c,"na12878",k,"all",n[k]}' >> $OUT/counts/$c.tsv
  echo "$(date -u +%FT%TZ) $c done"
}
export -f one_chrom; export OUT PANEL SAMPLES

CHROMS=$(bcftools index -s $PANEL | cut -f1)
/usr/bin/time -v xargs -P $SLURM_CPUS_PER_TASK -I{} bash -c 'one_chrom {}' <<< "$CHROMS"

{ printf 'chrom\tset\tclass\ttype\tmaf_bin\tn\n'; cat $OUT/counts/*.tsv; } > $OUT/counts-by-chrom.tsv
{ printf 'set\tclass\ttype\tmaf_bin\tn\n'
  awk -F'\t' 'NR>1{k=$2"\t"$3"\t"$4"\t"$5; n[k]+=$6} END{for(k in n) print k"\t"n[k]}' $OUT/counts-by-chrom.tsv | sort; } > $OUT/summary.tsv
column -t $OUT/summary.tsv
echo "end $(date -u +%FT%TZ)"
