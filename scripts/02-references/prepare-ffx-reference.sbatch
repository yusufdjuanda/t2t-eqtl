#!/bin/bash
#SBATCH --job-name=t2t-ffxref
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/storage/projects/yusuf.maulana/t2t-eqtl/logs/02-references/%x_%j.out
#SBATCH --error=/storage/projects/yusuf.maulana/t2t-eqtl/logs/02-references/%x_%j.err
# Stage 02 (references): prepare the third reference, FixItFelix's patched GRCh38 (Behera et al. 2023), for RNA-seq.
#
#   sbatch scripts/02-references/prepare-ffx-reference.sbatch
#
# Source: Zenodo record 7532965, GRCh38_masked_v2_decoy_gene.fasta.gz (md5 35d00c3c48d1849adb8e02086c3a1c70):
# the GRCh38 no-alt analysis set with the false copies masked (N) and decoy contigs for collapsed genes.
# Gene models: the same GRCh38 RefSeq RS_2025_08 GTF as grch38 (masked genes then receive no reads; decoys are
# unannotated, so reads placed there are not counted). Checks, written to compare-to-grch38.tsv:
#   every GRCh38 contig present with the same length; extra contigs (decoys) listed; N count per contig
#   (bases newly masked); chrM, chrEBV and masked chrY PARs as in grch38.
# Output: data/star-references/grch38-ffx/ in the same layout as the other references, with STAR indexes for
# 76 bp (star-2.7.11b-oh75, Geuvadis) and 150 bp reads (star-2.7.11b-oh149, MAGE).
set -euo pipefail
ROOT=/storage/home/yusuf.maulana/research/t2t-eqtl
URL="https://zenodo.org/records/7532965/files/GRCh38_masked_v2_decoy_gene.fasta.gz?download=1"
MD5=35d00c3c48d1849adb8e02086c3a1c70
G38=$ROOT/data/star-references/grch38
OUT=$ROOT/data/star-references/grch38-ffx
SRC=$ROOT/data/external/fixitfelix-behera2023
module load samtools/1.23.1 star/2.7.11b
export TMPDIR=/scratch/$USER/t2t-eqtl/tmp/$SLURM_JOB_ID; mkdir -p $TMPDIR; trap 'rm -rf $TMPDIR' EXIT
mkdir -p $OUT $SRC
echo "start $(date -u +%FT%TZ)"

F=$SRC/GRCh38_masked_v2_decoy_gene.fasta.gz
if [[ ! -f $F || -w $F ]]; then
  for i in $(seq 1 50); do curl -sSfL -C - --retry 5 --retry-all-errors -o $F "$URL" && break; sleep 20; done
  [[ $(md5sum < $F | cut -d' ' -f1) == $MD5 ]] || { echo "md5 mismatch for $F" >&2; exit 1; }
  chmod a-w $F
fi
echo "downloaded and md5-verified: $F"

zcat $F > $OUT/genome.fa
samtools faidx $OUT/genome.fa
ln -sfn $G38/annotation.gtf $OUT/annotation.gtf
cp $G38/genes.tsv $G38/tx2gene.tsv $G38/keep-regions.bed $OUT/

# compare with grch38: contig names and lengths, N counts per contig
pixi run --manifest-path $ROOT/pixi.toml python -I - $G38/genome.fa $OUT/genome.fa $OUT/compare-to-grch38.tsv <<'EOF'
import sys
def scan(fa):
    n, L, c = {}, {}, None
    for line in open(fa):
        if line[0] == '>':
            c = line[1:].split()[0]; n[c] = 0; L[c] = 0
        else:
            s = line.strip(); L[c] += len(s); n[c] += s.count('N') + s.count('n')
    return L, n
La, na = scan(sys.argv[1]); Lb, nb = scan(sys.argv[2])
with open(sys.argv[3], 'w') as o:
    o.write("contig\tin_grch38\tin_ffx\tlength_grch38\tlength_ffx\tN_grch38\tN_ffx\tnewly_masked\n")
    for c in list(La) + [c for c in Lb if c not in La]:
        o.write(f"{c}\t{c in La}\t{c in Lb}\t{La.get(c,'')}\t{Lb.get(c,'')}\t{na.get(c,'')}\t{nb.get(c,'')}\t"
                f"{nb.get(c,0)-na.get(c,0) if c in La and c in Lb else ''}\n")
missing = [c for c in La if c not in Lb]; lendiff = [c for c in La if c in Lb and La[c] != Lb[c]]
extra = [c for c in Lb if c not in La]
print(f"grch38 contigs {len(La)}, ffx {len(Lb)}; missing in ffx {missing}; length differs {lendiff}; extra (decoys) {extra}")
print("newly masked bases:", sum(nb[c]-na[c] for c in La if c in Lb), "on", [c for c in La if c in Lb and nb[c] != na[c]])
if missing or lendiff: sys.exit("ffx reference does not keep GRCh38 coordinates")
EOF

awk -F'\t' '$3!="gene"' $OUT/annotation.gtf > $TMPDIR/no-gene.gtf
pixi run --manifest-path $ROOT/pixi.toml gffread -w $OUT/transcripts.fa -g $OUT/genome.fa $TMPDIR/no-gene.gtf

for OH in 75 149; do
  IDX=$OUT/star-2.7.11b-oh$OH; rm -rf $IDX
  /usr/bin/time -v STAR --runMode genomeGenerate --runThreadN $SLURM_CPUS_PER_TASK --genomeDir $IDX \
    --genomeFastaFiles $OUT/genome.fa --sjdbGTFfile $OUT/annotation.gtf --sjdbOverhang $OH \
    --limitGenomeGenerateRAM 60000000000 --outTmpDir $TMPDIR/star-$OH --outFileNamePrefix $IDX/
done
( cd $OUT && md5sum genome.fa genome.fa.fai transcripts.fa genes.tsv tx2gene.tsv keep-regions.bed compare-to-grch38.tsv star-2.7.11b-oh*/* > CHECKSUMS.md5 )
du -sb $OUT
echo "end $(date -u +%FT%TZ)"
