"""Step 5: combine per-chromosome tensorQTL results and compute genome-wide q-values (Storey, R qvalue).

    pixi run python scripts/05-eqtl/collect-eqtl.py --cohort mage --ref chm13v2.0 [--method star] [--variants matched]

Reads results/05-eqtl/<cohort>/<ref>/<method>/<variants>/by-chrom/*.cis.tsv.gz; writes cis-eqtl.tsv.gz there with
tensorQTL's calculate_qvalues (q-values on the beta-approximated permutation p-values, eGenes at FDR 5%, and the
nominal p-value threshold per gene for use with the nominal pairs). Prints the number of eGenes.
"""
import argparse
from pathlib import Path

import pandas as pd
from tensorqtl import post

ROOT = Path(__file__).resolve().parents[2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--method", default="star")
    ap.add_argument("--variants", default="matched")
    ap.add_argument("--fdr", type=float, default=0.05)
    args = ap.parse_args()
    d = ROOT / f"results/05-eqtl/{args.cohort}/{args.ref}/{args.method}/{args.variants}"
    parts = sorted((d / "by-chrom").glob("*.cis.tsv.gz"))
    res = pd.concat(pd.read_csv(p, sep="\t", index_col=0) for p in parts)
    post.calculate_qvalues(res, fdr=args.fdr, qvalue_lambda=None)
    res.to_csv(d / "cis-eqtl.tsv.gz", sep="\t", float_format="%.6g")
    print(f"{args.cohort} {args.ref} {args.method} {args.variants}: {len(parts)} chromosomes, {len(res)} genes, "
          f"{int((res.qval <= args.fdr).sum())} eGenes at FDR {args.fdr}")


if __name__ == "__main__":
    main()
