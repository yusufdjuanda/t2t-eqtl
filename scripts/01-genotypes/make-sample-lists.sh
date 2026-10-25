#!/usr/bin/env bash
# Step 1: list which Geuvadis individuals have genotypes in the phased 1KGP T2T-CHM13v2.0 panel.
#
#   bash scripts/01-genotypes/make-sample-lists.sh
#
# Geuvadis individuals = unique sample_alias (GEUV:<id>, prefix removed) in the ENA run report for PRJEB3366
# (data/metadata/geuvadis/ena-PRJEB3366-2026-10-13.tsv; 667 runs, 464 individuals).
# Panel individuals = sample names in the panel BCF (3,202).
# Writes config/samples/geuvadis-in-1kgp-chm13.txt (both; 451) and geuvadis-not-in-1kgp-chm13.txt (Geuvadis only; 13).
# Light (a few seconds): reads only the BCF header. The lists were first made by hand on 2026-10-13; this script
# reproduces them exactly (checked 2026-10-17).
set -euo pipefail
ROOT=/storage/home/yusuf.maulana/research/t2t-eqtl
ENA=$ROOT/data/metadata/geuvadis/ena-PRJEB3366-2026-10-13.tsv
PANEL=$ROOT/data/external/1kgp-chm13v2.0-phased-shapeit5-v1.1/1KGP.CHM13v2.0.whole_genome.recalibrated.snp_indel.pass.phased.native_maps.3202.bcf.gz
OUT=$ROOT/config/samples
module load bcftools/1.23.1
tmp=$(mktemp -d -p /scratch/$USER/t2t-eqtl/tmp); trap 'rm -rf $tmp' EXIT

awk -F'\t' 'NR>1 {sub(/^GEUV:/, "", $2); print $2}' $ENA | sort -u > $tmp/geuvadis
bcftools query -l $PANEL | sort -u > $tmp/panel
comm -12 $tmp/geuvadis $tmp/panel > $OUT/geuvadis-in-1kgp-chm13.txt
comm -23 $tmp/geuvadis $tmp/panel > $OUT/geuvadis-not-in-1kgp-chm13.txt
echo "Geuvadis $(wc -l < $tmp/geuvadis), panel $(wc -l < $tmp/panel), in both $(wc -l < $OUT/geuvadis-in-1kgp-chm13.txt), Geuvadis only $(wc -l < $OUT/geuvadis-not-in-1kgp-chm13.txt)"
