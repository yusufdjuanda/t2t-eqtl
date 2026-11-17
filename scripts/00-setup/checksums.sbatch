#!/bin/bash
#SBATCH --job-name=t2t-sha256
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/00-setup/%x_%j.out
#SBATCH --error=logs/00-setup/%x_%j.err
# sha256 sums of the downloaded inputs as a Slurm job (see checksums.sh).
#
#   scripts/submit.sh scripts/00-setup/checksums.sbatch write|check
set -euo pipefail
MODE=${1:?write or check}
ROOT=${T2T_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && git rev-parse --show-toplevel)}}
source "$ROOT/scripts/lib/common.sh"
echo "start $(date -u +%FT%TZ) $MODE"
/usr/bin/time -v bash "$ROOT/scripts/00-setup/checksums.sh" "$MODE"
echo "end $(date -u +%FT%TZ)"
