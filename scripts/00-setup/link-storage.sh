#!/usr/bin/env bash
# Create the project's large folders (data, results, logs, work), either as links to project storage or as plain
# folders, and link the inputs shared with the sibling project t2t-atlas.
#
#   pixi run setup        (same as: bash scripts/00-setup/link-storage.sh)
#
# Reads config/site.env (template: config/site.env.example):
#   T2T_DATA_DIR, T2T_RESULTS_DIR, T2T_LOGS_DIR, T2T_WORK_DIR   if set, <root>/<folder> becomes a link to it;
#                                                                if unset, a plain folder (git ignores all four)
#   T2T_ATLAS_DIR   if set, data/references/<name> links to $T2T_ATLAS_DIR/data/references/<name> and
#                   data/external/t2t-atlas-analysis to $T2T_ATLAS_DIR/results/analysis (read-only inputs)
# Safe to rerun: existing correct links are left alone; a real folder is never replaced.
set -euo pipefail
ROOT=${T2T_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && git rev-parse --show-toplevel)}
source "$ROOT/scripts/lib/common.sh"

link() {  # link <target> <link path>
  local target=$1 path=$2
  if [[ -L $path ]]; then
    [[ $(readlink "$path") == "$target" ]] && { echo "ok      $path -> $target"; return; }
    echo "relink  $path -> $target (was $(readlink "$path"))"; rm "$path"
  elif [[ -e $path ]]; then
    echo "keep    $path (a real folder; not replaced by a link to $target)"; return
  fi
  [[ -e $target ]] || echo "warning: $target does not exist yet" >&2
  mkdir -p "$(dirname "$path")"; ln -s "$target" "$path"; echo "linked  $path -> $target"
}

for name in data results logs work; do
  var=T2T_${name^^}_DIR
  if [[ -n ${!var:-} ]]; then link "${!var}" "$ROOT/$name"
  else mkdir -p "$ROOT/$name"; echo "folder  $ROOT/$name"; fi
done

if [[ -n ${T2T_ATLAS_DIR:-} ]]; then
  for ref in chm13v2.0-ncbi grch38-noalt liftover t2t-annotations; do
    link "$T2T_ATLAS_DIR/data/references/$ref" "$ROOT/data/references/$ref"
  done
  link "$T2T_ATLAS_DIR/results/analysis" "$ROOT/data/external/t2t-atlas-analysis"
else
  echo "T2T_ATLAS_DIR not set: put the references in data/references/ yourself (see README, 'Input data')."
fi
mkdir -p "$ROOT/data/metadata" "$ROOT/data/external"
