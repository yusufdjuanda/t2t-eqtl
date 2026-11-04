#!/bin/bash
#SBATCH --job-name=t2t-staridx
#SBATCH --partition=cpu
#SBATCH --array=0-1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/storage/projects/yusuf.maulana/t2t-eqtl/logs/02-references/%x_%A_%a.out
#SBATCH --error=/storage/projects/yusuf.maulana/t2t-eqtl/logs/02-references/%x_%A_%a.err
# Stage 02 (references): build a STAR index for a given read length from an existing STAR reference (array 0 = chm13v2.0, 1 = grch38).
#
#   sbatch scripts/02-references/build-star-index.sbatch <sjdbOverhang>      e.g. 149 for MAGE (2 x 150 bp)
#
# Reuses genome.fa and annotation.gtf from data/star-references/<ref>/ (made by
# scripts/02-references/prepare-star-references.sbatch); only the splice-junction overhang differs.
# Output: data/star-references/<ref>/star-2.7.11b-oh<overhang>/ (the pilot index star-2.7.11b-oh75/ is overhang 75).
set -euo pipefail
OVERHANG=${1:?usage: build-star-index.sbatch <sjdbOverhang>}
ROOT=/storage/home/yusuf.maulana/research/t2t-eqtl
IDS=(chm13v2.0 grch38)
ID=${IDS[$SLURM_ARRAY_TASK_ID]}
REF=$ROOT/data/star-references/$ID
IDX=$REF/star-2.7.11b-oh$OVERHANG
module load star/2.7.11b
export TMPDIR=/scratch/$USER/t2t-eqtl/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}; mkdir -p $TMPDIR; trap 'rm -rf $TMPDIR' EXIT
echo "start $(date -u +%FT%TZ) $ID overhang $OVERHANG"
rm -rf $IDX
/usr/bin/time -v STAR --runMode genomeGenerate --runThreadN $SLURM_CPUS_PER_TASK \
  --genomeDir $IDX --genomeFastaFiles $REF/genome.fa --sjdbGTFfile $REF/annotation.gtf --sjdbOverhang $OVERHANG \
  --limitGenomeGenerateRAM 60000000000 --outTmpDir $TMPDIR/star --outFileNamePrefix $IDX/
( cd $REF && md5sum star-2.7.11b-oh$OVERHANG/* > CHECKSUMS-star-2.7.11b-oh$OVERHANG.md5 )
du -sb $IDX
echo "end $(date -u +%FT%TZ) $ID"
