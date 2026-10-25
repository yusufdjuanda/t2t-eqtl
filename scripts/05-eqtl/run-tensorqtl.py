"""Step 5: cis-eQTL mapping with tensorQTL for one chromosome (called by run-tensorqtl.sbatch).

    pixi run python scripts/05-eqtl/run-tensorqtl.py --cohort mage --ref chm13v2.0 --method star --variants matched \
        --chrom chr22 [--nperm 10000] [--threads 16]

Inputs: results/05-eqtl/<cohort>/<ref>/<method>/{phenotypes.bed.gz, covariates.tsv} and the genotypes of the reference
(chm13v2.0 -> genotypes/t2t.<variants>; grch38 and grch38-ffx -> genotypes/grch38.<variants>).
Runs, for the genes on this chromosome (TSS +- 1 Mb, in-cohort MAF >= 0.01):
  map_cis       permutations with beta approximation -> one top variant per gene with empirical p-values
  map_nominal   all variant-gene pairs (parquet)
Writes results/05-eqtl/<cohort>/<ref>/<method>/<variants>/by-chrom/<chrom>.cis.tsv.gz and nominal/<chrom>.*.parquet.
q-values across chromosomes are computed by collect-eqtl.py.
"""
import argparse
from pathlib import Path

import pandas as pd
import torch
from tensorqtl import cis, pgen

ROOT = Path(__file__).resolve().parents[2]
GENO = {"chm13v2.0": "t2t", "grch38": "grch38", "grch38-ffx": "grch38"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--method", default="star")
    ap.add_argument("--variants", default="matched", choices=["matched", "all"])
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--nperm", type=int, default=10000)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()
    torch.set_num_threads(args.threads)

    d = ROOT / f"results/05-eqtl/{args.cohort}/{args.ref}/{args.method}"
    out = d / args.variants
    (out / "by-chrom").mkdir(parents=True, exist_ok=True)
    (out / "nominal").mkdir(parents=True, exist_ok=True)

    ph = pd.read_csv(d / "phenotypes.bed.gz", sep="\t")
    ph = ph[ph["#chr"] == args.chrom]
    if ph.empty:
        print(f"{args.chrom}: no genes"); return
    pos = ph[["#chr", "end", "phenotype_id"]].rename(columns={"#chr": "chr", "end": "pos"}).set_index("phenotype_id")
    ph = ph.set_index("phenotype_id").iloc[:, 3:]
    cov = pd.read_csv(d / "covariates.tsv", sep="\t", index_col=0).T

    reader = pgen.PgenReader(str(ROOT / f"results/05-eqtl/{args.cohort}/genotypes/{GENO[args.ref]}.{args.variants}"))
    variant_df = reader.variant_df
    variant_df = variant_df[variant_df.chrom == args.chrom]
    geno = reader.read_list(variant_df.index.tolist(), impute_mean=True)
    geno = pd.DataFrame(geno, index=variant_df.index, columns=reader.sample_ids)
    common = [s for s in ph.columns if s in geno.columns and s in cov.index]
    ph, geno, cov = ph[common], geno[common], cov.loc[common]
    print(f"{args.cohort} {args.ref} {args.method} {args.variants} {args.chrom}: {ph.shape[0]} genes, "
          f"{geno.shape[0]} variants, {len(common)} individuals, {cov.shape[1]} covariates")

    res = cis.map_cis(geno, variant_df[["chrom", "pos"]], ph, pos, covariates_df=cov, maf_threshold=0.01,
                      nperm=args.nperm, seed=args.seed, verbose=False)
    res.to_csv(out / "by-chrom" / f"{args.chrom}.cis.tsv.gz", sep="\t", float_format="%.6g")
    cis.map_nominal(geno, variant_df[["chrom", "pos"]], ph, pos, f"{args.chrom}", covariates_df=cov,
                    maf_threshold=0.01, output_dir=str(out / "nominal"), verbose=False)


if __name__ == "__main__":
    main()
