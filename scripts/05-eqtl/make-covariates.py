"""Step 5: covariates for tensorQTL: genotype PCs (shared by all references), sex, expression PCs (per reference).

    pixi run python scripts/05-eqtl/make-covariates.py --cohort mage --ref chm13v2.0 --method star [--n-expr-pcs auto]

Genotype PCs are computed ONCE per cohort (on first use) from the T2T *matched* variants (results/05-eqtl/<cohort>/
genotypes/t2t.matched.*): MAF >= 0.05, LD-pruned (plink2 --indep-pairwise 200kb 1 0.2), plink2 --pca 5. Because only
matched variants are used, the same individuals' ancestry is described identically for every reference.
Sex: from config/samplesheets/<cohort>.csv (male = 1, female = 0).
Expression PCs: PCA of the phenotype BED of this reference and method. Number: --n-expr-pcs, or 'auto' = the GTEx rule
by sample size (15 for < 150 individuals, 30 for 150-249, 45 for 250-349, 60 for >= 350); FixItFelix chose 30 for
449 individuals by eGene saturation, which can be checked later with --n-expr-pcs.
Writes results/05-eqtl/<cohort>/<ref>/<method>/covariates.tsv (covariates x individuals, tensorQTL format).
"""
import argparse
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def genotype_pcs(cohort):
    g = ROOT / f"results/05-eqtl/{cohort}/genotypes"
    out = g / "genotype-pcs.eigenvec"
    if not out.exists():
        run = lambda *a: subprocess.run(["pixi", "run", "--manifest-path", str(ROOT / "pixi.toml"), "plink2", *a],
                                        check=True, stdout=subprocess.DEVNULL)
        run("--pfile", str(g / "t2t.matched"), "--maf", "0.05", "--snps-only", "--autosome",
            "--indep-pairwise", "200kb", "1", "0.2", "--out", str(g / "genotype-prune"))
        run("--pfile", str(g / "t2t.matched"), "--extract", str(g / "genotype-prune.prune.in"),
            "--pca", "5", "--out", str(g / "genotype-pcs"))
    pcs = pd.read_csv(out, sep="\t").rename(columns={"#IID": "IID"}).set_index("IID")
    return pcs[[c for c in pcs.columns if c.startswith("PC")]].T.rename(index=lambda c: f"geno_{c}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--method", default="star")
    ap.add_argument("--n-expr-pcs", default="auto")
    args = ap.parse_args()
    d = ROOT / f"results/05-eqtl/{args.cohort}/{args.ref}/{args.method}"
    ph = pd.read_csv(d / "phenotypes.bed.gz", sep="\t", index_col=3).iloc[:, 3:]
    n = ph.shape[1]
    k = int(args.n_expr_pcs) if args.n_expr_pcs != "auto" else (15 if n < 150 else 30 if n < 250 else 45 if n < 350 else 60)
    k = min(k, n - 10)
    x = ph.T.values - ph.T.values.mean(axis=0)
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    epc = pd.DataFrame((u[:, :k] * s[:k]).T, columns=ph.columns, index=[f"expr_PC{i + 1}" for i in range(k)])
    sheet = pd.read_csv(ROOT / f"config/samplesheets/{args.cohort}.csv").drop_duplicates("sample").set_index("sample")
    sex = pd.DataFrame([sheet.sex.reindex(ph.columns).map({"male": 1, "female": 0}).values], columns=ph.columns, index=["sex"])
    gpc = genotype_pcs(args.cohort).reindex(columns=ph.columns)
    cov = pd.concat([gpc, sex, epc])
    missing = cov.columns[cov.isna().any()]
    if len(missing):
        print(f"warning: {len(missing)} individuals lack a covariate (sex or genotype PCs), filled with the mean: {list(missing)[:5]}")
        cov = cov.apply(lambda r: r.fillna(r.mean()), axis=1)
    cov.index.name = "covariate"
    cov.to_csv(d / "covariates.tsv", sep="\t", float_format="%.6g")
    var = s[:k] ** 2 / (s ** 2).sum()
    print(f"{args.cohort} {args.ref} {args.method}: {n} individuals; {gpc.shape[0]} genotype PCs, sex, {k} expression PCs "
          f"({100 * var.sum():.1f}% of expression variance)")


if __name__ == "__main__":
    main()
