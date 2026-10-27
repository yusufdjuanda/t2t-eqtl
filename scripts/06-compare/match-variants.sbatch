#!/bin/bash
#SBATCH --job-name=t2t-varmatch
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=23
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=/storage/projects/yusuf.maulana/t2t-eqtl/logs/06-compare/%x_%j.out
#SBATCH --error=/storage/projects/yusuf.maulana/t2t-eqtl/logs/06-compare/%x_%j.err
# Step 6: pair common variants between the T2T and GRCh38 genotype sets of the 1,089 expression individuals.
#
#   sbatch scripts/06-compare/match-variants.sbatch [min MAF, default 0.01]
#
# Inputs: results/01-genotypes/chm13v2.0/expression/<chr>.bcf (T2T) and
# results/01-genotypes/grch38/expression/<chr>.bcf (GRCh38). Both: multiallelic sites split,
# MAF recomputed in these individuals. Liftover (UCSC hs1ToHg38 chain) is used ONLY to pair variants that were called
# independently on each reference, never to create genotypes.
# Per chromosome (chr1-22, X), variants with MAF >= min in either set:
#   1. T2T variant -> 1-bp interval -> liftOver to GRCh38 (keeps the chain strand)
#   2. match on GRCh38 chrom + position + alleles: same REF/ALT ("same"), or REF/ALT swapped because the two
#      references carry different bases ("swapped"); on a minus-strand chain alleles are reverse-complemented first
#   3. categories: t2t: matched (MAF within 0.02) | matched_maf_discordant | no_liftover (T2T-only sequence) |
#                  lifted_not_called_grch38 ;  grch38: not_called_t2t (GRCh38 variants with no T2T partner)
# Output: results/06-compare/variant-matching/<chr>.tsv.gz (all variants, both sides, category, MAF on each side,
# T2T SYNTENIC flag) and summary.tsv (counts by category x type x MAF bin).
set -euo pipefail
MINMAF=${1:-0.01}
ROOT=/storage/home/yusuf.maulana/research/t2t-eqtl
T2T=$ROOT/results/01-genotypes/chm13v2.0/expression
G38=$ROOT/results/01-genotypes/grch38/expression
CHAIN=$ROOT/data/references/liftover/hs1ToHg38.over.chain.gz
OUT=$ROOT/results/06-compare/variant-matching
module load bcftools/1.23.1
export TMPDIR=/scratch/$USER/t2t-eqtl/tmp/$SLURM_JOB_ID; mkdir -p $TMPDIR; trap 'rm -rf $TMPDIR' EXIT
mkdir -p $OUT
echo "start $(date -u +%FT%TZ) min MAF $MINMAF"

one_chrom() {
  local c=$1 d=$TMPDIR/$1; mkdir -p $d
  bcftools query -i "INFO/MAF>=$MINMAF" -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/MAF\t%INFO/SYNTENIC\n' $T2T/$c.bcf > $d/t2t.tsv
  bcftools query -i "INFO/MAF>=$MINMAF" -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/MAF\n' $G38/$c.bcf > $d/g38.tsv
  awk -F'\t' -v OFS='\t' '{print $1, $2-1, $2, NR, 0, "+"}' $d/t2t.tsv > $d/t2t.bed
  pixi run --manifest-path $ROOT/pixi.toml liftOver -minMatch=0.95 $d/t2t.bed $CHAIN $d/lifted.bed $d/unlifted.bed 2>/dev/null
  pixi run --manifest-path $ROOT/pixi.toml python -I $ROOT/scripts/06-compare/match-variants.py $d $OUT/$c.tsv.gz
  echo "$(date -u +%FT%TZ) $c done"
}
export -f one_chrom; export T2T G38 CHAIN OUT ROOT MINMAF TMPDIR

/usr/bin/time -v xargs -P $SLURM_CPUS_PER_TASK -I{} bash -c 'set -euo pipefail; one_chrom {}' < <(printf 'chr%s\n' {1..22} X)
pixi run --manifest-path $ROOT/pixi.toml python -I - $OUT <<'EOF'
import sys, glob, pandas as pd
out = sys.argv[1]
d = pd.concat(pd.read_csv(f, sep="\t") for f in sorted(glob.glob(f"{out}/chr*.tsv.gz")))
d["maf"] = d[["maf_t2t", "maf_grch38"]].max(axis=1)
d["maf_bin"] = pd.cut(d.maf, [0, 0.05, 0.5001], labels=["0.01-0.05", ">=0.05"], right=False)
s = d.groupby(["side", "category", "type", "maf_bin"], observed=True).size().rename("n").reset_index()
s.to_csv(f"{out}/summary.tsv", sep="\t", index=False)
print(s.pivot_table(index=["side", "category"], columns=["type", "maf_bin"], values="n", aggfunc="sum", observed=True).fillna(0).astype(int).to_string())
EOF
echo "end $(date -u +%FT%TZ)"
