#!/usr/bin/env bash
# Step 1: which individuals of the expression cohorts (MAGE, Geuvadis) have genotypes in the phased 1KGP
# T2T-CHM13v2.0 panel, i.e. can be used for eQTL mapping.
#
#   bash scripts/01-genotypes/make-sample-lists.sh
#
# Sources (all in data/metadata/ and data/external/; how they were obtained: docs/t2t-eqtl.md 4.2-4.3):
#   MAGE      ENA run report for PRJNA851328 (mage/ena-PRJNA851328-2026-10-17.tsv): sample_alias = individual
#   Geuvadis  ENA run report for PRJEB3366 (geuvadis/ena-PRJEB3366-fastq-md5-2026-10-17.tsv): sample_alias without
#             the GEUV: prefix = individual
#   panel     sample names in the panel BCF header (3,202; bcftools query -l reads the header only)
# Writes config/samples/expression-individuals.tsv, one row per individual in either cohort:
#   sample  mage  geuvadis  t2t_panel  included      (yes/no; included = in a cohort and in the panel)
# Read by scripts/03-expression/make-cohort-sheets.py (which individuals' runs go into the samplesheets) and by the
# genotype subset scripts (01-genotypes/chm13v2.0/subset-genotypes.sbatch, 01-genotypes/grch38/subset-genotypes.sbatch).

set -euo pipefail
ROOT=/storage/home/yusuf.maulana/research/t2t-eqtl
MAGE=$ROOT/data/metadata/mage/ena-PRJNA851328-2026-10-17.tsv
GEUV=$ROOT/data/metadata/geuvadis/ena-PRJEB3366-fastq-md5-2026-10-17.tsv
PANEL=$ROOT/data/external/1kgp-chm13v2.0-phased-shapeit5-v1.1/1KGP.CHM13v2.0.whole_genome.recalibrated.snp_indel.pass.phased.native_maps.3202.bcf.gz
OUT=$ROOT/config/samples/expression-individuals.tsv
module load bcftools/1.23.1
mkdir -p /scratch/$USER/t2t-eqtl/tmp
tmp=$(mktemp -d -p /scratch/$USER/t2t-eqtl/tmp); trap 'rm -rf $tmp' EXIT

# individual IDs from the sample_alias column, found by name
alias_col() { awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="sample_alias") c=i; next} {sub(/^GEUV:/, "", $c); print $c}' "$1" | sort -u; }
alias_col $MAGE > $tmp/mage
alias_col $GEUV > $tmp/geuvadis
bcftools query -l $PANEL | sort -u > $tmp/panel

sort -u $tmp/mage $tmp/geuvadis \
  | awk -v m=$tmp/mage -v g=$tmp/geuvadis -v p=$tmp/panel 'BEGIN{
        OFS = "\t"
        while ((getline x < m) > 0) M[x] = 1
        while ((getline x < g) > 0) G[x] = 1
        while ((getline x < p) > 0) P[x] = 1
        print "sample", "mage", "geuvadis", "t2t_panel", "included" }
      function yn(v) { return v ? "yes" : "no" }
      { print $1, yn($1 in M), yn($1 in G), yn($1 in P), yn($1 in P) }' > $OUT

awk -F'\t' 'NR>1{
    if ($2=="yes") { m++; if ($4=="yes") mi++ }
    if ($3=="yes") { g++; if ($4=="yes") gi++ }
    if ($2=="yes" && $3=="yes") b++
    if ($5=="yes") n++ }
  END{ printf "MAGE %d (%d in panel), Geuvadis %d (%d in panel), in both cohorts %d; included %d\n", m, mi, g, gi, b, n }' $OUT
