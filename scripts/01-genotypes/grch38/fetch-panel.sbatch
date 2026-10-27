#!/bin/bash
#SBATCH --job-name=t2t-g38geno
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=2-00:00:00
#SBATCH --output=/storage/projects/yusuf.maulana/t2t-eqtl/logs/01-genotypes/grch38/%x_%j.out
#SBATCH --error=/storage/projects/yusuf.maulana/t2t-eqtl/logs/01-genotypes/grch38/%x_%j.err
# Step 4: download the phased 1000 Genomes GRCh38 panel (NYGC high coverage, 3,202 samples; Byrska-Bishop et al.
# 2022): SNVs, indels and SVs, chr1-22 and chrX (v2), with tabix indexes.
#
#   sbatch scripts/01-genotypes/grch38/fetch-panel.sbatch
#
# Source: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/
#         20220422_3202_phased_SNV_INDEL_SV/
# Plain HTTP, because HTTPS to this server and to the AWS mirror is reset from Karakoram (2026-10-17). Integrity
# comes from the md5s in the server's 20220804_manifest.txt: every file is checked, and a verified file is made
# read-only and skipped on rerun. Downloads resume (-C -) and retry. 4 parallel streams.
set -euo pipefail
ROOT=/storage/home/yusuf.maulana/research/t2t-eqtl
URL=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV
OUT=$ROOT/data/external/1kgp-grch38-nygc-phased-20220422
export TMPDIR=/scratch/$USER/t2t-eqtl/tmp/$SLURM_JOB_ID; mkdir -p $TMPDIR; trap 'rm -rf $TMPDIR' EXIT
mkdir -p $OUT; cd $OUT
echo "start $(date -u +%FT%TZ)"
for f in 20220804_manifest.txt change_log.txt README_1kGP_phased_panel_110722.pdf; do
  [[ -s $f ]] || curl -sSf --retry 10 --retry-all-errors -o $f $URL/$f
done

fetch_one() {  # file md5
  local f=$1 md5=$2
  if [[ -f $f && ! -w $f ]]; then echo "$f already verified"; return 0; fi
  for i in $(seq 1 100); do
    curl -sSf -C - --retry 5 --retry-all-errors --retry-delay 20 -o $f $URL/$f && break
    echo "$(date -u +%T) $f attempt $i failed at $(stat -c %s $f 2>/dev/null || echo 0) bytes" >&2; sleep 30
  done
  if [[ $(md5sum < $f | cut -d' ' -f1) == $md5 ]]; then chmod a-w $f; echo "$(date -u +%T) $f md5 ok"
  else echo "$f md5 MISMATCH" >&2; return 1; fi
}
export -f fetch_one; export URL

/usr/bin/time -v xargs -P 4 -L1 bash -c 'set -euo pipefail; fetch_one "$@"' _ < 20220804_manifest.txt
du -sb $OUT
echo "end $(date -u +%FT%TZ)"
