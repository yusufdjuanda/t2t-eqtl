#!/usr/bin/env bash
# Static checks run by CI (.github/workflows/check.yml) and by `pixi run -e lint lint`:
#   1. every Python script compiles (syntax);
#   2. every shell and Slurm script parses (bash -n);
#   3. ShellCheck finds no errors (warnings and style notes are not enforced);
#   4. no script hard-codes a cluster path (paths come from scripts/lib/common.sh and config/site.env);
#   5. the data-sources table is well formed.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
fail=0

echo "1. Python syntax"
python -m compileall -q scripts || fail=1

echo "2. shell syntax"
mapfile -t sh < <(git ls-files 'scripts/*.sh' 'scripts/*.sbatch')
for f in "${sh[@]}"; do bash -n "$f" || { echo "   $f"; fail=1; }; done

echo "3. ShellCheck (errors only)"
shellcheck --shell=bash --severity=error "${sh[@]}" || fail=1

echo "4. hard-coded cluster paths"
if git grep -nE '/storage/|/scratch/|/home/[a-z]' -- 'scripts/*' 'config/*.env.example' ':!scripts/00-setup/lint.sh'; then
  echo "   use \$ROOT, \$T2T_TMP or \$T2T_SCRATCH instead"; fail=1
fi

echo "5. config/data-sources.tsv"
awk -F'\t' 'NR == 1 {n = NF} NF != n {print "   line " NR ": " NF " columns, expected " n; bad = 1} END {exit bad}' \
  config/data-sources.tsv || fail=1

[[ $fail -eq 0 ]] && echo "all checks passed" || { echo "some checks failed"; exit 1; }
