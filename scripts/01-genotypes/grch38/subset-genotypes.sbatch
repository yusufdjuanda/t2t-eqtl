#!/bin/bash
#SBATCH --job-name=t2t-g38sub
#SBATCH --cpus-per-task=23
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/01-genotypes/grch38/%x_%j.out
#SBATCH --error=logs/01-genotypes/grch38/%x_%j.err
# Step 4: subset the phased 1000 Genomes GRCh38 panel (NYGC) to the expression individuals, matching the T2T subset
# (scripts/01-genotypes/chm13v2.0/subset-genotypes.sbatch) as closely as the panels allow.
#
#   scripts/submit.sh scripts/01-genotypes/grch38/subset-genotypes.sbatch [output folder] [sample list]
#
# Same filters as the T2T subset: the given individuals; multiallelic sites split (norm -m-); AC/AN/MAF recomputed
# in this subset; monomorphic sites dropped. One addition: SNVs and indels only (-v snps,indels), because the NYGC
# panel also holds SVs and the T2T panel does not. Per chromosome (chr1-22, chrX), in parallel:
#   expression/<chr>.bcf   the individuals in the sample list (default: the included individuals of
#                          config/samples/expression-individuals.tsv, 1,089 MAGE + Geuvadis)
#   na12878/<chr>.bcf      NA12878, sites where it carries a non-reference allele
#   counts-by-chrom.tsv    variant counts by type x MAF bin, and NA12878 het / hom-alt counts
# Output: summary.tsv (all chromosomes), default folder results/01-genotypes/grch38/.
set -euo pipefail
ROOT=${T2T_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && git rev-parse --show-toplevel)}}
source "$ROOT/scripts/lib/common.sh"
OUT=${1:-$ROOT/results/01-genotypes/grch38}
SAMPLES=${2:-}
PANEL_DIR=$ROOT/data/external/1kgp-grch38-nygc-phased-20220422
t2t_env bio   # STAR, salmon, samtools, ... from pixi.toml (published runs: cluster modules of the same versions)
t2t_tmpdir "$SLURM_JOB_ID"
mkdir -p $OUT/expression $OUT/na12878 $TMPDIR/counts
# default sample list: included individuals of config/samples/expression-individuals.tsv (make-sample-lists.sh)
if [ -z "$SAMPLES" ]; then
  SAMPLES=$TMPDIR/samples.txt
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="included") c=i; next} $c=="yes"{print $1}' \
    $ROOT/config/samples/expression-individuals.tsv > $SAMPLES
fi
echo "start $(date -u +%FT%TZ)"

one_chrom() {
  local c=$1 f
  f=$(ls $PANEL_DIR/1kGP_high_coverage_Illumina.$c.filtered.SNV_INDEL_SV_phased_panel*.vcf.gz)
  bcftools view -S $SAMPLES -v snps,indels -Ou $f 2>/dev/null \
    | bcftools norm -m- -Ou \
    | bcftools +fill-tags -Ou -- -t AC,AN,MAF \
    | bcftools view -c1:minor -Ob -W=csi -o $OUT/expression/$c.bcf
  bcftools view -s NA12878 -v snps,indels -c1 -Ob -W=csi -o $OUT/na12878/$c.bcf $f 2>/dev/null
  bcftools query -f '%TYPE\t%INFO/MAF\n' $OUT/expression/$c.bcf \
    | awk -v c=$c 'BEGIN{OFS="\t"} {t=($1=="SNP")?"snp":"indel";
        b=($2>=0.05)?"maf>=0.05":($2>=0.01)?"maf0.01-0.05":"maf<0.01"; n[t"\t"b]++}
        END{for(k in n) print c,"expression",k,n[k]}' > $TMPDIR/counts/$c.tsv
  bcftools query -f '[%GT]\n' $OUT/na12878/$c.bcf \
    | awk -v c=$c 'BEGIN{OFS="\t"} {split($1,g,/[|\/]/); z=(g[1]!=g[2])?"het":"hom_alt"; n[z]++}
        END{for(k in n) print c,"na12878",k,"all",n[k]}' >> $TMPDIR/counts/$c.tsv
  echo "$(date -u +%FT%TZ) $c done"
}
export -f one_chrom; export OUT SAMPLES PANEL_DIR TMPDIR

CHROMS=$(printf 'chr%s\n' {1..22} X)
/usr/bin/time -v xargs -P $SLURM_CPUS_PER_TASK -I{} bash -c 'set -euo pipefail; one_chrom {}' <<< "$CHROMS"

{ printf 'chrom\tset\ttype\tmaf_bin\tn\n'; cat $TMPDIR/counts/*.tsv; } > $OUT/counts-by-chrom.tsv
{ printf 'set\ttype\tmaf_bin\tn\n'
  awk -F'\t' 'NR>1{k=$2"\t"$3"\t"$4; n[k]+=$5} END{for(k in n) print k"\t"n[k]}' $OUT/counts-by-chrom.tsv | sort; } > $OUT/summary.tsv
column -t $OUT/summary.tsv
echo "end $(date -u +%FT%TZ)"
