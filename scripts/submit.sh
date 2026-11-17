#!/usr/bin/env bash
# Submit a Slurm batch script of this project from the project root, with the site's partition and account.
#
#   scripts/submit.sh [sbatch options] <scripts/.../job.sbatch> [job arguments]
#   scripts/submit.sh --dependency=afterok:123 scripts/05-eqtl/run-tensorqtl.sbatch mage chm13v2.0
#   pixi run submit scripts/06-compare/coloc.sbatch mage
#
# Why: batch scripts write their logs to logs/<stage>/ (relative paths in their #SBATCH lines) and find the project
# through SLURM_SUBMIT_DIR, so they must be submitted from the root. Partitions and accounts differ between clusters,
# so they are not written in the scripts: a script that needs the long partition says so in a comment line
# "# t2t-partition: long", and config/site.env maps "default" and "long" to this cluster's partition names
# (T2T_PARTITION, T2T_PARTITION_LONG) and sets T2T_ACCOUNT if needed. Options given here override both.
# Prints sbatch's output (with --parsable: the job ID only).
set -euo pipefail
ROOT=${T2T_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && git rev-parse --show-toplevel)}
source "$ROOT/scripts/lib/common.sh"
cd "$ROOT"

opts=()
while [[ $# -gt 0 && $1 == -* ]]; do opts+=("$1"); shift; done
job=${1:?usage: scripts/submit.sh [sbatch options] <job.sbatch> [arguments]}; shift
[[ -f $job ]] || { echo "error: $job not found (paths are relative to the project root)" >&2; exit 1; }

want=$(sed -n 's/^# t2t-partition: *\([a-z]*\).*/\1/p' "$job" | head -1)
case ${want:-default} in
  long) part=${T2T_PARTITION_LONG:-${T2T_PARTITION:-}} ;;
  *)    part=${T2T_PARTITION:-} ;;
esac
site=()
[[ -n $part ]] && site+=(--partition="$part")
[[ -n ${T2T_ACCOUNT:-} ]] && site+=(--account="$T2T_ACCOUNT")

# Slurm does not create the log folder; make it from the script's --output line.
logdir=$(sed -n 's/^#SBATCH --output=\(.*\)\/[^/]*$/\1/p' "$job" | head -1)
[[ -n $logdir ]] && mkdir -p "$logdir"

exec sbatch "${site[@]}" "${opts[@]}" "$job" "$@"
