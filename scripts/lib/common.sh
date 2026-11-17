# shellcheck shell=bash
# Shared set-up for every shell and Slurm script in this project. Source it right after locating the project root:
#
#   ROOT=${T2T_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && git rev-parse --show-toplevel)}}
#   source "$ROOT/scripts/lib/common.sh"
#
# Why the root is found this way: Slurm copies a batch script to a spool folder before running it, so the script's own
# path does not point into the project. Jobs are therefore submitted from the project root (scripts/submit.sh does
# this), and SLURM_SUBMIT_DIR is the root. Plain shell scripts find the root through git. T2T_ROOT overrides both.
#
# After sourcing:
#   T2T_ROOT      project root (absolute)
#   T2T_SCRATCH   large temporary space (default: <root>/work); set in config/site.env
#   T2T_TMP       parent of every job's TMPDIR (default: $T2T_SCRATCH/tmp)
#   t2t_env [env] puts a pixi environment's tools on PATH (default env: default). This replaces `module load`.
#   t2t_tmpdir <name>  makes and exports TMPDIR=$T2T_TMP/<name> and removes it when the script exits.

if [[ -z ${ROOT:-} || ! -f $ROOT/pixi.toml || ! -f $ROOT/scripts/lib/common.sh ]]; then
  echo "error: project root not found (ROOT='${ROOT:-}'). Submit jobs from the project root with" \
       "scripts/submit.sh, or set T2T_ROOT." >&2
  exit 1
fi
ROOT=$(cd "$ROOT" && pwd -P)
export T2T_ROOT=$ROOT

# Site settings (scratch location, Slurm partitions). config/site.env is not in git; copy config/site.env.example.
if [[ -f $ROOT/config/site.env ]]; then
  set -a
  # shellcheck source=/dev/null
  source "$ROOT/config/site.env"
  set +a
fi
export T2T_SCRATCH=${T2T_SCRATCH:-$ROOT/work}
export T2T_TMP=${T2T_TMP:-$T2T_SCRATCH/tmp}

# Put the tools of a pixi environment (pixi.toml) first on PATH in the current shell.
# Only PATH changes, on purpose: a full `pixi shell-hook` activation would make every later `pixi run python ...` in the
# same script use this environment instead of the default one.
t2t_env() {
  local env=${1:-default}
  if [[ ! -x $ROOT/.pixi/envs/$env/bin/python && ! -d $ROOT/.pixi/envs/$env/conda-meta ]]; then
    pixi install --locked --manifest-path "$ROOT/pixi.toml" --environment "$env" >&2 || {
      echo "error: pixi environment '$env' is not installed (run: pixi install --locked -e $env)" >&2; exit 1; }
  fi
  export PATH=$ROOT/.pixi/envs/$env/bin:$PATH
}

# A private temporary folder for this job, removed on exit. TMPDIR is used by sort, samtools, STAR, Python, R, ...
# Never use /tmp on the cluster: it is the node's small shared root disk.
t2t_tmpdir() {
  export TMPDIR=$T2T_TMP/${1:?t2t_tmpdir needs a name}
  mkdir -p "$TMPDIR"
  # shellcheck disable=SC2064
  trap "rm -rf '$TMPDIR'" EXIT
}
