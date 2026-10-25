#!/usr/bin/env bash
# Step 6: for every gene, list the other genes whose exons overlap its exons, on T2T and on GRCh38, and flag genes
# whose overlap partners differ between the references.
#
#   bash scripts/06-compare/compare-annotation-overlaps.sh      (a few minutes; bedtools on the two GTFs)
#
# Why: STAR GeneCounts discards a read as ambiguous when it overlaps exons of more than one gene (Geuvadis and MAGE are
# unstranded, so strand does not separate them). If a gene has an exon-overlapping neighbour on one reference but not
# the other, its STAR count changes between references for annotation reasons alone, while salmon (which assigns
# such reads to transcripts) changes much less. Pilots: about half of the changed genes in class `other` changed in
# STAR only (docs/t2t-eqtl.md 13A.4b).
# Output: results/06-compare/annotation-overlaps.tsv, one row per GeneID annotated on both references:
#   gene_ncbi  partners_t2t  partners_grch38  gained_on_t2t  lost_on_t2t  overlap_differs
# (partners = GeneIDs of exon-overlapping genes, ';'-separated; grch38-ffx uses the grch38 annotation).
set -euo pipefail
ROOT=/storage/home/yusuf.maulana/research/t2t-eqtl
OUT=$ROOT/results/06-compare
module load bedtools/2.31.1
mkdir -p $OUT /scratch/$USER/t2t-eqtl/tmp
TMP=$(mktemp -d -p /scratch/$USER/t2t-eqtl/tmp); trap 'rm -rf "$TMP"' EXIT

for ref in chm13v2.0 grch38; do
  # exons -> BED (chrom, start, end, gene_id), sorted
  awk -F'\t' -v OFS='\t' '$3=="exon" { g=$9; sub(/.*gene_id "/, "", g); sub(/".*/, "", g); print $1, $4-1, $5, g }' \
      $ROOT/data/star-references/$ref/annotation.gtf \
    | sort -k1,1 -k2,2n -S 4G -T $TMP > $TMP/$ref.exons.bed
  # gene_id -> GeneID
  awk -F'\t' 'NR>1 {print $1"\t"$2}' $ROOT/data/star-references/$ref/genes.tsv > $TMP/$ref.ids
  bedtools intersect -a $TMP/$ref.exons.bed -b $TMP/$ref.exons.bed -wa -wb -sorted \
    | awk -F'\t' '$4 != $8 {print $4"\t"$8}' | sort -u -S 4G -T $TMP > $TMP/$ref.pairs
  awk -F'\t' 'NR==FNR {id[$1]=$2; next} ($1 in id) && ($2 in id) && id[$1] != id[$2] {print id[$1]"\t"id[$2]}' \
    $TMP/$ref.ids $TMP/$ref.pairs | sort -u > $TMP/$ref.idpairs
  echo "$ref: $(wc -l < $TMP/$ref.exons.bed) exons, $(wc -l < $TMP/$ref.idpairs) overlapping GeneID pairs"
done

pixi run --manifest-path $ROOT/pixi.toml python -I - $TMP $ROOT <<'EOF'
import sys
import pandas as pd
tmp, root = sys.argv[1], sys.argv[2]
def partners(ref):
    p = pd.read_csv(f"{tmp}/{ref}.idpairs", sep="\t", header=None, names=["a", "b"], dtype=str)
    return p.groupby("a").b.apply(set)
t, g = partners("chm13v2.0"), partners("grch38")
ids = lambda ref: set(pd.read_csv(f"{root}/data/star-references/{ref}/genes.tsv", sep="\t", dtype={"gene_ncbi": str}).gene_ncbi)
both = sorted(ids("chm13v2.0") & ids("grch38"))
rows = []
for i in both:
    a, b = t.get(i, set()), g.get(i, set())
    rows.append((i, ";".join(sorted(a)), ";".join(sorted(b)), ";".join(sorted(a - b)), ";".join(sorted(b - a)), a != b))
d = pd.DataFrame(rows, columns=["gene_ncbi", "partners_t2t", "partners_grch38", "gained_on_t2t", "lost_on_t2t", "overlap_differs"])
d.to_csv(f"{root}/results/06-compare/annotation-overlaps.tsv", sep="\t", index=False)
print(f"genes on both references: {len(d)}; exon-overlap partners differ: {d.overlap_differs.sum()} "
      f"(gained on T2T {(d.gained_on_t2t != '').sum()}, lost on T2T {(d.lost_on_t2t != '').sum()})")
EOF
