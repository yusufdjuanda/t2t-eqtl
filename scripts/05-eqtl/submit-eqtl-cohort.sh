#!/usr/bin/env bash
# Steps 3-6 for one cohort after quantification: submit the whole eQTL chain as Slurm jobs with dependencies.
#
#   bash scripts/05-eqtl/submit-eqtl-cohort.sh <cohort> [methods=star,salmon] [variants=matched] [chromosomes=chr1..chr22]
#
#   merge-expression.sbatch                         matrices (step 3)
#   prepare-genotypes.sbatch  t2t, grch38            pgen sets (after merge: individuals = matrix columns)
#   prepare-phenotypes-covariates.sbatch per method  (after both genotype jobs; methods one after another)
#   run-tensorqtl.sbatch per reference and method    arrays over the chromosomes, at most 2 tasks each at a time
#   collect-eqtl.sbatch per reference and method
#   compare-eqtl.sbatch per method                   (after the three collects of that method)
#   combine-methods.py                               (after all compares, when both star and salmon are run)
# Steps with afterok dependencies: a failure stops what depends on it (scontrol show job <id> shows DependencyNeverSatisfied).
# Job IDs go to logs/05-eqtl/submitted-<cohort>-<time>.tsv.
set -euo pipefail
COHORT=${1:?cohort}; METHODS=${2:-star,salmon}; VARIANTS=${3:-matched}
CHROMS=${4:-$(printf 'chr%s,' {1..22} | sed 's/,$//')}
ROOT=/storage/home/yusuf.maulana/research/t2t-eqtl
NCHR=$(tr ',' '\n' <<< "$CHROMS" | wc -l)
[[ $CHROMS == $(printf 'chr%s,' $(seq 1 $NCHR) | sed 's/,$//') ]] || { echo "chromosomes must be chr1..chrN (array index = chromosome - 1)" >&2; exit 1; }
LOG=$ROOT/logs/05-eqtl/submitted-$COHORT-$(date -u +%Y%m%dT%H%M%SZ).tsv
printf 'step\targs\tjob\n' > $LOG
sub() { local step=$1; shift; local j; j=$(sbatch --parsable "$@"); printf '%s\t%s\t%s\n' "$step" "${*: -4}" "$j" >> $LOG; echo $j; }

M=$(sub merge $ROOT/scripts/04-expression/merge-expression.sbatch $COHORT)
G1=$(sub genotypes --dependency=afterok:$M $ROOT/scripts/05-eqtl/prepare-genotypes.sbatch $COHORT t2t $CHROMS)
G2=$(sub genotypes --dependency=afterok:$M $ROOT/scripts/05-eqtl/prepare-genotypes.sbatch $COHORT grch38 $CHROMS)
CMP=(); PREVP=""
for METHOD in ${METHODS//,/ }; do
  # one at a time: the first computes the cohort's genotype PCs, the others reuse them
  P=$(sub phenocov --dependency=afterok:$G1:$G2${PREVP:+:$PREVP} $ROOT/scripts/05-eqtl/prepare-phenotypes-covariates.sbatch $COHORT $METHOD); PREVP=$P
  C=()
  for REF in chm13v2.0 grch38 grch38-ffx; do
    T=$(sub tensorqtl --dependency=afterok:$P --array=0-$((NCHR-1))%2 $ROOT/scripts/05-eqtl/run-tensorqtl.sbatch $COHORT $REF $METHOD $VARIANTS)
    C+=($(sub collect --dependency=afterok:$T $ROOT/scripts/05-eqtl/collect-eqtl.sbatch $COHORT $REF $METHOD $VARIANTS))
  done
  CMP+=($(sub compare --dependency=afterok:$(IFS=:; echo "${C[*]}") $ROOT/scripts/06-compare/compare-eqtl.sbatch $COHORT $METHOD $VARIANTS))
done
if [[ $METHODS == *star* && $METHODS == *salmon* ]]; then
  sub combine --dependency=afterok:$(IFS=:; echo "${CMP[*]}") --job-name=t2t-combine --partition=cpu --cpus-per-task=2 --mem=16G \
    --time=01:00:00 --output=/storage/projects/yusuf.maulana/t2t-eqtl/logs/06-compare/%x_%j.out \
    --error=/storage/projects/yusuf.maulana/t2t-eqtl/logs/06-compare/%x_%j.err \
    --wrap="bash -c 'cd $ROOT && pixi run python scripts/06-compare/combine-methods.py --cohort $COHORT --variants $VARIANTS'" >/dev/null
fi
cat $LOG
