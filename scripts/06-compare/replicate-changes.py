"""Step 6: replication of the robust changed genes of a discovery cohort in a replication cohort, at the same variants.

    pixi run python scripts/06-compare/replicate-changes.py [--discovery mage] [--replication geuvadis] [--variants matched]

Discovery set: results/06-compare/eqtl-comparison/<discovery>/robust-<variants>/changed-<A>-vs-<B>.tsv without artefact
flag (the 57-gene set for MAGE, T2T vs patched GRCh38). Each gene has a STAR lead variant pair (var_a on reference A,
var_b on B, the same variant in both panels) and an allele alignment sign_b = sign(geno_r).

Why the same variants: z_diff = (slope_A - slope_B) / se is signed by the lead variant's allele coding, and in about
half of the genes the slopes have opposite signs on the two references, so "stronger_on" (larger |slope|) and the
sign of z_diff from separately chosen leads are not comparable between cohorts. Here the replication cohort's nominal
statistics (map_nominal parquet, STAR and salmon) are read for exactly var_a and var_b, aligned with the discovery
sign_b, and z_rep is compared with the discovery STAR z_diff; the variant IDs and alleles are the same in both cohorts.
  tested        both variants present (in-cohort MAF >= 0.01) for the gene under both methods
  same_sign     sign(z_rep) == sign(z_disc) under STAR and salmon
  nominal       same_sign and p < 0.05 under both methods
  bonferroni    same_sign and p < 0.05 / (genes tested) under both methods
A sign test (binomial, one-sided) on the STAR z_rep asks whether same sign exceeds 50%.
Writes results/06-compare/eqtl-comparison/replication-<replication>/{genes-<A>-vs-<B>.tsv, summary.tsv}.

--exclude-shared (2026-10-28): 93 of Geuvadis's 451 individuals are also in MAGE (same genotypes, separate RNA-seq),
so the replication is not fully independent. With this flag the individuals of the discovery cohort are dropped and
the nominal statistics are recomputed on the rest (ols_nominal, the tensorQTL map_nominal computation on the step-5
phenotypes and covariates, which were normalised and computed on all individuals); outputs go to
replication-<replication>-independent/. --check compares ols_nominal on all individuals with the parquet values.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from tensorqtl import pgen

ROOT = Path(__file__).resolve().parents[2]
PAIRS = [("chm13v2.0", "grch38"), ("chm13v2.0", "grch38-ffx"), ("grch38-ffx", "grch38")]
METHODS = ("star", "salmon")
GENO = {"chm13v2.0": "t2t", "grch38": "grch38", "grch38-ffx": "grch38"}
_cache = {}


def ols_nominal(cohort, ref, method, variants, want, exclude=()):
    """Nominal slope, slope_se, pval_nominal for (phenotype_id, variant_id) pairs in want, recomputed without the
    individuals in exclude, as tensorQTL map_nominal: genotype (missing = mean) and phenotype residualised on the
    covariates plus intercept, dof = n - 2 - covariates."""
    key = (cohort, ref, method, variants)
    if key not in _cache:
        d = ROOT / f"results/05-eqtl/{cohort}/{ref}/{method}"
        ph = pd.read_csv(d / "phenotypes.bed.gz", sep="\t").set_index("phenotype_id").iloc[:, 3:]
        cov = pd.read_csv(d / "covariates.tsv", sep="\t", index_col=0).T
        rd = pgen.PgenReader(str(ROOT / f"results/05-eqtl/{cohort}/genotypes/{GENO[ref]}.{variants}"))
        _cache[key] = ph, cov, rd
    ph, cov, rd = _cache[key]
    samples = [s for s in rd.sample_ids if s in cov.index and s in ph.columns and s not in set(exclude)]
    sidx = [rd.sample_ids.index(s) for s in samples]
    C = np.column_stack([np.ones(len(samples)), cov.loc[samples].values])
    dof = len(samples) - C.shape[1] - 1
    Q, _ = np.linalg.qr(C)
    have = set(rd.variant_ids)
    rows = []
    for r in want.itertuples():
        if r.phenotype_id not in ph.index or r.variant_id not in have:
            continue
        g = np.asarray(rd.read(r.variant_id), dtype=float)[sidx]
        g = np.where(np.isnan(g), np.nanmean(g), g)
        y = ph.loc[r.phenotype_id, samples].values.astype(float)
        g, y = g - Q @ (Q.T @ g), y - Q @ (Q.T @ y)
        slope = (g @ y) / (g @ g)
        se = np.sqrt(((y - slope * g) @ (y - slope * g)) / dof / (g @ g))
        rows.append((r.phenotype_id, r.variant_id, slope, se, 2 * stats.t.sf(abs(slope / se), dof)))
    out = pd.DataFrame(rows, columns=["phenotype_id", "variant_id", "slope", "slope_se", "pval_nominal"])
    out.attrs["n"] = len(samples)
    return out.set_index(["phenotype_id", "variant_id"])


def nominal(cohort, ref, method, variants, want):
    """Nominal stats for (phenotype_id, variant_id) pairs in want (DataFrame with chrom, phenotype_id, variant_id)."""
    d = ROOT / f"results/05-eqtl/{cohort}/{ref}/{method}/{variants}/nominal"
    out = []
    for chrom, w in want.groupby("chrom"):
        f = d / f"{chrom}.cis_qtl_pairs.{chrom}.parquet"
        n = pd.read_parquet(f, columns=["phenotype_id", "variant_id", "slope", "slope_se"],
                            filters=[("phenotype_id", "in", w.phenotype_id.unique().tolist())])
        out.append(n.merge(w[["phenotype_id", "variant_id"]], on=["phenotype_id", "variant_id"]))
    return pd.concat(out).set_index(["phenotype_id", "variant_id"])


def shared_individuals(discovery, replication):
    """Individuals with RNA-seq in both cohorts (runs.tsv sample IDs)."""
    s = [set(pd.read_csv(ROOT / f"results/04-expression/{c}/matrices/runs.tsv", sep="\t")["sample"])
         for c in (discovery, replication)]
    return s[0] & s[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--discovery", default="mage")
    ap.add_argument("--replication", default="geuvadis")
    ap.add_argument("--variants", default="matched")
    ap.add_argument("--exclude-shared", action="store_true", help="drop the discovery cohort's individuals")
    ap.add_argument("--check", action="store_true", help="recompute on all individuals and compare with tensorQTL")
    args = ap.parse_args()
    base = ROOT / "results/06-compare/eqtl-comparison"
    out = base / f"replication-{args.replication}{'-independent' if args.exclude_shared else ''}"
    exclude = shared_individuals(args.discovery, args.replication) if args.exclude_shared else set()
    if args.exclude_shared:
        print(f"excluding {len(exclude)} {args.replication} individuals also in {args.discovery}")
    if not args.check:
        out.mkdir(parents=True, exist_ok=True)
    gid = {r: pd.read_csv(ROOT / f"data/star-references/{r}/genes.tsv", sep="\t", dtype={"gene_ncbi": str})
           .drop_duplicates("gene_ncbi").set_index("gene_ncbi").gene_id for r in ("chm13v2.0", "grch38", "grch38-ffx")}
    rows = []
    for a, b in PAIRS:
        t = pd.read_csv(base / f"{args.discovery}/robust-{args.variants}/changed-{a}-vs-{b}.tsv", sep="\t",
                        dtype={"gene_ncbi": str}).set_index("gene_ncbi")
        t = t[~t.artefact_flag & t.var_a.notna() & t.var_b.notna()]
        t = t[["symbol", "region_class", "var_a", "var_b", "geno_r", "slope_a", "slope_b", "z_diff"]].rename(
            columns={"slope_a": "disc_slope_a", "slope_b": "disc_slope_b", "z_diff": "disc_z"})
        t["sign_b"] = np.sign(t.geno_r)
        for m in METHODS:
            st = {}
            for ref, var in ((a, "var_a"), (b, "var_b")):
                w = pd.DataFrame({"phenotype_id": gid[ref].reindex(t.index).values, "variant_id": t[var].values,
                                  "chrom": t[var].str.split(":").str[0].values}).dropna()
                if args.check:
                    x = nominal(args.replication, ref, m, args.variants, w)
                    y = ols_nominal(args.replication, ref, m, args.variants, w).reindex(x.index)
                    print(f"check {ref} {m}: {len(x)} pairs, max |slope diff| {(x.slope - y.slope).abs().max():.2g}, "
                          f"max relative se diff {((x.slope_se - y.slope_se) / x.slope_se).abs().max():.2g}")
                    continue
                st[ref] = (ols_nominal(args.replication, ref, m, args.variants, w, exclude) if args.exclude_shared
                           else nominal(args.replication, ref, m, args.variants, w))
            if args.check:
                continue
            ka = pd.MultiIndex.from_arrays([gid[a].reindex(t.index), t.var_a])
            kb = pd.MultiIndex.from_arrays([gid[b].reindex(t.index), t.var_b])
            sa, sb = st[a].reindex(ka), st[b].reindex(kb)
            t[f"{m}_slope_a"] = sa.slope.values
            t[f"{m}_slope_b"] = sb.slope.values * t.sign_b.values
            t[f"{m}_z"] = (t[f"{m}_slope_a"] - t[f"{m}_slope_b"]) / np.sqrt(sa.slope_se.values ** 2 + sb.slope_se.values ** 2)
            t[f"{m}_p"] = 2 * stats.norm.sf(t[f"{m}_z"].abs())
        if args.check:
            continue
        t["tested"] = t.star_z.notna() & t.salmon_z.notna()
        n = int(t.tested.sum())
        t["same_sign"] = t.tested & (np.sign(t.star_z) == np.sign(t.disc_z)) & (np.sign(t.salmon_z) == np.sign(t.disc_z))
        pmax = t[["star_p", "salmon_p"]].max(axis=1)
        t["nominal"] = t.same_sign & (pmax < 0.05)
        t["bonferroni"] = t.same_sign & (pmax < 0.05 / max(n, 1))
        t.sort_values(["bonferroni", "nominal", "star_p"], ascending=[False, False, True]).to_csv(
            out / f"genes-{a}-vs-{b}.tsv", sep="\t", float_format="%.4g")
        s = t[t.star_z.notna()]
        k = int((np.sign(s.star_z) == np.sign(s.disc_z)).sum())
        rows.append({"pair": f"{a} vs {b}", "discovery_genes": len(t), "tested": n, "same_sign": int(t.same_sign.sum()),
                     "nominal": int(t.nominal.sum()), "bonferroni": int(t.bonferroni.sum()),
                     "star_tested": len(s), "star_same_sign": k,
                     "sign_test_p": stats.binomtest(k, len(s), 0.5, alternative="greater").pvalue if len(s) else np.nan})
        print(f"\n{a} vs {b}: {len(t)} discovery genes, {n} tested at the same variants in {args.replication}; not replicated:")
        print(t.loc[t.tested & ~t.nominal, ["symbol", "region_class", "disc_z", "star_z", "star_p", "salmon_z",
                                            "salmon_p"]].round(4).to_string())
    if args.check:
        return
    s = pd.DataFrame(rows)
    s.to_csv(out / "summary.tsv", sep="\t", index=False, float_format="%.3g")
    pd.set_option("display.width", 250)
    print("\n" + s.to_string(index=False))


if __name__ == "__main__":
    main()
