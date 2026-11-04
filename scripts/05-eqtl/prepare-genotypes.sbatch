#!/bin/bash
#SBATCH --job-name=t2t-genoprep
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/storage/projects/yusuf.maulana/t2t-eqtl/logs/05-eqtl/%x_%j.out
#SBATCH --error=/storage/projects/yusuf.maulana/t2t-eqtl/logs/05-eqtl/%x_%j.err
# Step 5: genotypes for tensorQTL (plink2 pgen) for one cohort, from the T2T or the GRCh38 genotype subset.
#
#   sbatch scripts/05-eqtl/prepare-genotypes.sbatch <cohort> <t2t|grch38> [chromosomes, default chr1..chr22,chrX]
#
# Individuals: the columns of the cohort's expression matrices (results/04-expression/<cohort>/matrices/chm13v2.0/
# star-counts.tsv.gz), i.e. one per person. Variants: in-cohort MAF >= 0.01; IDs chr:pos:ref:alt (long alleles
# truncated), the same IDs as the variant-matching tables. Two sets:
#   <geno>.all.{pgen,pvar,psam}       all variants passing MAF
#   <geno>.matched.{pgen,pvar,psam}   only variants paired between references (results/06-compare/variant-matching/),
#                                     the set for the main comparison (docs/t2t-eqtl.md 13B)
# Inputs: t2t    results/01-genotypes/chm13v2.0/expression/<chr>.bcf
#         grch38 results/01-genotypes/grch38/expression/<chr>.bcf
# Output: results/05-eqtl/<cohort>/genotypes/.  The grch38 genotypes serve both grch38 and grch38-ffx.
set -euo pipefail
COHORT=${1:?cohort}; GENO=${2:?t2t|grch38}
CHROMS=${3:-$(printf 'chr%s,' {1..22} X | sed 's/,$//')}
ROOT=/storage/home/yusuf.maulana/research/t2t-eqtl
case $GENO in
  t2t)    SRC=$ROOT/results/01-genotypes/chm13v2.0/expression; IDCOL=t2t_id ;;
  grch38) SRC=$ROOT/results/01-genotypes/grch38/expression; IDCOL=grch38_id ;;
  *) echo "geno must be t2t or grch38" >&2; exit 1 ;;
esac
OUT=$ROOT/results/05-eqtl/$COHORT/genotypes
MATCH=$ROOT/results/06-compare/variant-matching
module load bcftools/1.23.1
export TMPDIR=/scratch/$USER/t2t-eqtl/tmp/$SLURM_JOB_ID; mkdir -p $TMPDIR; trap 'rm -rf $TMPDIR' EXIT
mkdir -p $OUT
echo "start $(date -u +%FT%TZ) $COHORT $GENO $CHROMS"

# first line only, without a pipe that would fail under pipefail when head closes it early
python3 -c "import gzip,sys; print('\\n'.join(gzip.open(sys.argv[1],'rt').readline().rstrip('\\n').split('\\t')[1:]))" \
  $ROOT/results/04-expression/$COHORT/matrices/chm13v2.0/star-counts.tsv.gz > $OUT/individuals.txt
echo "$(wc -l < $OUT/individuals.txt) individuals"
# sex for plink2 (needed for chrX): samplesheet sex male/female -> 1/2, missing -> 0 (unknown)
awk -F, 'NR==1 {for (i=1;i<=NF;i++) c[$i]=i; print "#IID\tSEX"; next}
  !seen[$c["sample"]]++ {x=($c["sex"]=="male")?1:($c["sex"]=="female")?2:0; print $c["sample"]"\t"x}' \
  $ROOT/config/samplesheets/$COHORT.csv > $OUT/sex.txt
ls $(for c in ${CHROMS//,/ }; do echo $SRC/$c.bcf; done) > $TMPDIR/files.txt
bcftools concat -f $TMPDIR/files.txt --threads $SLURM_CPUS_PER_TASK -Ob -o $TMPDIR/all.bcf
/usr/bin/time -v pixi run --manifest-path $ROOT/pixi.toml plink2 --bcf $TMPDIR/all.bcf --keep $OUT/individuals.txt --update-sex $OUT/sex.txt \
  --maf 0.01 --output-chr chr26 --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 200 truncate --rm-dup force-first \
  --threads $SLURM_CPUS_PER_TASK --memory 48000 --make-pgen --out $OUT/$GENO.all

# matched variants: IDs of this side's matched variants in the variant-matching tables
pixi run --manifest-path $ROOT/pixi.toml python -I - $MATCH $IDCOL $TMPDIR/matched.ids "$CHROMS" <<'EOF'
import sys, pandas as pd
mdir, col, out, chroms = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4].split(",")
ids = []
for c in chroms:
    d = pd.read_csv(f"{mdir}/{c}.tsv.gz", sep="\t", usecols=["side", "category", "t2t_id", "grch38_id"], low_memory=False)
    ids += d.loc[(d.side == "t2t") & (d.category == "matched"), col].dropna().tolist()
pd.Series(ids).to_csv(out, index=False, header=False)
print(len(ids), "matched variant IDs")
EOF
pixi run --manifest-path $ROOT/pixi.toml plink2 --pfile $OUT/$GENO.all --extract $TMPDIR/matched.ids --output-chr chr26 \
  --threads $SLURM_CPUS_PER_TASK --make-pgen --out $OUT/$GENO.matched
for s in all matched; do echo "$GENO.$s: $(grep -vc '^#' $OUT/$GENO.$s.pvar) variants"; done
echo "end $(date -u +%FT%TZ)"
