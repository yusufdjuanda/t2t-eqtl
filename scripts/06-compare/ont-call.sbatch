#!/bin/bash
#SBATCH --job-name=t2t-ontcall
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/storage/projects/yusuf.maulana/t2t-eqtl/logs/06-compare/%x_%A_%a.out
#SBATCH --error=/storage/projects/yusuf.maulana/t2t-eqtl/logs/06-compare/%x_%A_%a.err
# Step 6: small variants from 1KGP ONT long reads in the ASE target regions, one individual per array task.
#
#   sbatch --array=1-<n>%20 scripts/06-compare/ont-call.sbatch
#
# Individuals: the ONT individuals of results/06-compare/ase/targets/individuals.tsv (ont = True), task i = i-th.
# Per individual: stream the reads in targets/ont-t2t.bed from the remote T2T CRAM (1KG_ONT_VIENNA/t2t/<id>.t2t.cram;
# only the needed slices are fetched; retried for the connection resets seen from Karakoram) with the consortium's
# T2T reference (data/external/1kg-ont-vienna/1KG_ONT_VIENNA_t2t.fa) -> Clair3 (model r941_prom_sup_g5014: R9.4.1,
# Guppy 6 sup, as the 1KGP ONT basecalling) in those regions, with read-backed phasing (WhatsHap, --enable_phasing).
# Output: results/06-compare/ase/ont/<id>/{reads.bam(.bai), clair3.vcf.gz(.tbi), clair3.phased.vcf.gz(.tbi)},
# complete.txt (individuals that have it are skipped).
set -euo pipefail
ROOT=/storage/home/yusuf.maulana/research/t2t-eqtl
T=$ROOT/results/06-compare/ase/targets
ID=$(awk -F'\t' 'NR > 1 && $2 == "True" {print $1}' $T/individuals.tsv | sed -n "${SLURM_ARRAY_TASK_ID}p")
[[ -n $ID ]] || { echo "no individual for task $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
OUT=$ROOT/results/06-compare/ase/ont/$ID
[[ -f $OUT/complete.txt ]] && { echo "$ID complete, skipped"; exit 0; }
REF=$ROOT/data/external/1kg-ont-vienna/1KG_ONT_VIENNA_t2t.fa
URL=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/t2t/$ID.t2t.cram
export TMPDIR=/scratch/$USER/t2t-eqtl/tmp/ont-$SLURM_ARRAY_JOB_ID-$ID; mkdir -p $TMPDIR; trap 'rm -rf $TMPDIR' EXIT
mkdir -p $OUT
cd $TMPDIR
echo "start $(date -u +%FT%TZ) $ID"
px() { pixi run --manifest-path $ROOT/pixi.toml -e longread "$@"; }

for i in $(seq 1 20); do
  rc=0; px samtools view -@ $SLURM_CPUS_PER_TASK -b -M -L $T/ont-t2t.bed -T $REF -o reads.unsorted.bam $URL 2> view.err || rc=$?
  [[ $rc -eq 0 ]] && break
  echo "$(date -u +%T) fetch attempt $i exit $rc: $(head -c 300 view.err)" >&2; sleep 30
done
[[ $rc -eq 0 ]] || { echo "$ID: fetch failed" >&2; exit 1; }
px samtools sort -@ $SLURM_CPUS_PER_TASK -o $OUT/reads.bam reads.unsorted.bam
px samtools index $OUT/reads.bam
echo "$(date -u +%T) reads: $(px samtools view -c -F 0x900 $OUT/reads.bam) primary"

px bash -c "run_clair3.sh --bam_fn=$OUT/reads.bam --ref_fn=$REF --threads=$SLURM_CPUS_PER_TASK --platform=ont \
  --model_path=\$CONDA_PREFIX/bin/models/r941_prom_sup_g5014 --bed_fn=$T/ont-t2t.bed --sample_name=$ID \
  --enable_phasing --output=$TMPDIR/clair3" > clair3.log 2>&1 || { tail -30 clair3.log >&2; exit 1; }
cp clair3/merge_output.vcf.gz $OUT/clair3.vcf.gz; cp clair3/merge_output.vcf.gz.tbi $OUT/clair3.vcf.gz.tbi
cp clair3/phased_merge_output.vcf.gz $OUT/clair3.phased.vcf.gz; cp clair3/phased_merge_output.vcf.gz.tbi $OUT/clair3.phased.vcf.gz.tbi
echo "$(date -u +%T) variants: $(px bcftools view -H -f PASS $OUT/clair3.vcf.gz | wc -l) PASS"
echo "end $(date -u +%FT%TZ)"
date -u +%FT%TZ > $OUT/complete.txt
