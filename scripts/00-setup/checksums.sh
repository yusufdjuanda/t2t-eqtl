#!/usr/bin/env bash
# sha256 sums of every downloaded input, so a rebuilt project can prove it used the same files.
#
#   pixi run checksums        write  config/checksums/inputs.sha256 (after downloading; ~60 GB, ~10 min on 8 CPUs)
#   pixi run verify-data      check  every file listed there (prints only failures, then a summary)
#   scripts/submit.sh scripts/00-setup/checksums.sbatch write|check     the same as a Slurm job
#
# Covered: data/external/ (genotype panels, ONT assemblies, FixItFelix, GWAS, ctyper, T2T consortium BED) and
# data/references/ (genomes, annotations, liftOver chains; links to t2t-atlas on the original cluster are followed).
# Not covered: data/external/t2t-atlas-analysis (another project's results), derived files (STAR indexes, results/).
# Paths in the list are relative to the project root. Sources and URLs: config/data-sources.tsv.
set -euo pipefail
ROOT=${T2T_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && git rev-parse --show-toplevel)}}
source "$ROOT/scripts/lib/common.sh"
cd "$ROOT"
MODE=${1:?usage: checksums.sh write|check}
LIST=config/checksums/inputs.sha256
JOBS=${SLURM_CPUS_PER_TASK:-$(nproc)}

case $MODE in
  write)
    mkdir -p "$(dirname "$LIST")"
    find -L data/external data/references -type f -not -path 'data/external/t2t-atlas-analysis/*' -print0 \
      | sort -z | xargs -0 -n 1 -P "$JOBS" sha256sum | sort -k 2 > "$LIST.tmp"
    mv "$LIST.tmp" "$LIST"
    echo "wrote $LIST: $(wc -l < "$LIST") files"
    ;;
  check)
    [[ -f $LIST ]] || { echo "error: $LIST not found (run: pixi run checksums)" >&2; exit 1; }
    # one sha256sum -c per file, in parallel; failures and missing files are printed
    fails=$(awk '{print $2}' "$LIST" | xargs -n 1 -P "$JOBS" -I{} sh -c \
      'grep -F "  {}" "$0" | sha256sum -c --quiet 2>&1 || true' "$LIST" | tee /dev/stderr | grep -c ': FAILED' || true)
    echo "checked $(wc -l < "$LIST") files: ${fails:-0} failed"
    [[ ${fails:-0} -eq 0 ]]
    ;;
  *) echo "usage: checksums.sh write|check" >&2; exit 1 ;;
esac
