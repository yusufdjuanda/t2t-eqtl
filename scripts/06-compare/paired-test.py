"""Step 6: paired test of the eQTL effect change between references, accounting for the shared individuals.

    pixi run python scripts/06-compare/paired-test.py --cohort mage [--variants matched]

Why: compare-effects.py uses z_diff = (slope_A - slope_B) / sqrt(se_A^2 + se_B^2), which treats the two slopes as
independent. They are estimated from the same individuals and (matched) genotypes, and most of a gene's expression is
measured the same way on both references, so the two estimates are strongly positively correlated and that z is
conservative. Here, per gene, at the variant pair compare-effects.py chose (var_a, var_b, sign_b = sign(geno_r)):
  slopes and se are recomputed as tensorQTL map_nominal does (genotype, missing = mean, and expression residualised on
  each reference's own covariates plus intercept; dof = n - 2 - covariates), and
    cov(slope_A, slope_B) = [sum(e_A e_B) / dof] * sum(g_A g_B) / (sum(g_A^2) sum(g_B^2))
  (e = residuals of each reference's regression, g = residualised genotypes), the covariance of two OLS slopes whose
  errors are correlated within an individual, so
    z_paired = (slope_A - slope_B) / sqrt(se_A^2 + se_B^2 - 2 cov).
  rho = corr(e_A, e_B) is reported: how much of the unexplained expression the two references share.
Calibration: for --nboot-genes random genes, the individuals are resampled --nboot times and the SD of slope_A -
slope_B is compared with the paired and the independent standard errors (calibration.tsv).
Calls per method: changed = Bonferroni p_paired < 0.05 over the genes compared. Robust (the project rule, as in
combine-methods.py) = changed under STAR and salmon with the same sign. The 57-gene set came from the independent z,
so the output reports how the paired test adds to it.
Inputs: results/06-compare/eqtl-comparison/<cohort>/<method>-<variants>/effects-<A>-vs-<B>.tsv and the step-5
phenotypes, covariates and genotypes. Writes results/06-compare/eqtl-comparison/<cohort>/paired-<variants>/
{effects-<A>-vs-<B>-<method>.tsv, robust-<A>-vs-<B>.tsv, calibration.tsv, summary.tsv}.

--copy-number: only genes with a ctyper copy number (results/06-compare/ctyper/copy-number.tsv.gz, the gene and the
family copy number where they vary), which are added as covariates for that gene on both references (projected out of
expression and genotype after the step-5 covariates; dof reduced accordingly). Adds cn_r2 (variance of the lead
genotype explained by the copy numbers, on reference A). Output in paired-<variants>-cn/.

--at <discovery>: replication. Instead of the cohort's own lead pairs, the unflagged robust genes of the discovery
cohort's paired run (paired-<variants>/robust-<A>-vs-<B>.tsv) are tested at the discovery variant pairs (same IDs and
alleles; sign from the discovery geno_r), without the individuals also in the discovery cohort. Replicated = same sign as the discovery STAR z_paired under STAR and salmon
with p < 0.05 (nominal) or p < 0.05 / genes tested (Bonferroni). Output in paired-<variants>-at-<discovery>/
{effects-*.tsv, replication-<A>-vs-<B>.tsv, replication-summary.tsv}.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm
from tensorqtl import pgen

ROOT = Path(__file__).resolve().parents[2]
PAIRS = [("chm13v2.0", "grch38"), ("chm13v2.0", "grch38-ffx"), ("grch38-ffx", "grch38")]
METHODS = ("star", "salmon")
GENO = {"chm13v2.0": "t2t", "grch38": "grch38", "grch38-ffx": "grch38"}


def resid(X, Q):
    """Rows of X residualised on the orthonormal basis Q (individuals x covariates)."""
    return X - (X @ Q) @ Q.T


def replication(res, out):
    """Replication calls for the discovery genes tested at the discovery variant pairs (--at)."""
    rows = []
    for a, b in PAIRS:
        s, q = res[("star", a, b)], res[("salmon", a, b)]
        j = s.join(q[["delta", "z_paired", "p_paired"]], rsuffix="_salmon", how="inner")
        n = len(j)
        same = (np.sign(j.z_paired) == np.sign(j.disc_z)) & (np.sign(j.z_paired_salmon) == np.sign(j.disc_z))
        pmax = j[["p_paired", "p_paired_salmon"]].max(axis=1)
        j["same_sign"], j["nominal"], j["bonferroni"] = same, same & (pmax < 0.05), same & (pmax < 0.05 / max(n, 1))
        j.sort_values("p_paired").to_csv(out / f"replication-{a}-vs-{b}.tsv", sep="\t", float_format="%.4g")
        big = j.disc_delta.abs() >= 0.1
        for name, k in (("all", j.index == j.index), ("in independent-z set", j.in_indep_set.astype(bool)),
                        ("new with paired test", ~j.in_indep_set.astype(bool)),
                        ("new, |delta| >= 0.1", ~j.in_indep_set.astype(bool) & big),
                        ("new, |delta| < 0.1", ~j.in_indep_set.astype(bool) & ~big)):
            x = j[k]
            rows.append({"pair": f"{a} vs {b}", "genes": name, "tested": len(x), "same_sign": int(x.same_sign.sum()),
                         "nominal": int(x.nominal.sum()), "bonferroni": int(x.bonferroni.sum())})
    s = pd.DataFrame(rows)
    s.to_csv(out / "replication-summary.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print("\n" + s.to_string(index=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--variants", default="matched")
    ap.add_argument("--nboot", type=int, default=200)
    ap.add_argument("--nboot-genes", type=int, default=300)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--copy-number", action="store_true", help="add ctyper copy numbers as per-gene covariates")
    ap.add_argument("--at", help="discovery cohort: test its paired robust genes at its variant pairs")
    args = ap.parse_args()
    base = ROOT / f"results/06-compare/eqtl-comparison/{args.cohort}"
    out = base / f"paired-{args.variants}{'-cn' if args.copy_number else ''}{f'-at-{args.at}' if args.at else ''}"
    cnt = None
    shared = set()
    if args.at:
        shared = [set(pd.read_csv(ROOT / f"results/04-expression/{c}/matrices/runs.tsv", sep="\t")["sample"])
                  for c in (args.at, args.cohort)]
        shared = shared[0] & shared[1]
        print(f"excluding {len(shared)} individuals also in {args.at}")
    if args.copy_number:
        cnt = pd.read_csv(ROOT / "results/06-compare/ctyper/copy-number.tsv.gz", sep="\t")
        cnt = cnt[cnt.iloc[:, 2:].var(axis=1) > 0]
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    gid = {r: pd.read_csv(ROOT / f"data/star-references/{r}/genes.tsv", sep="\t", dtype={"gene_ncbi": str})
           .drop_duplicates("gene_ncbi").set_index("gene_ncbi").gene_id for r in ("chm13v2.0", "grch38", "grch38-ffx")}
    rd = {s: pgen.PgenReader(str(ROOT / f"results/05-eqtl/{args.cohort}/genotypes/{s}.{args.variants}"))
          for s in ("t2t", "grch38")}

    res, calib = {}, []
    for m in METHODS:
        ph, cov = {}, {}
        for r in ("chm13v2.0", "grch38", "grch38-ffx"):
            d = ROOT / f"results/05-eqtl/{args.cohort}/{r}/{m}"
            ph[r] = pd.read_csv(d / "phenotypes.bed.gz", sep="\t").set_index("phenotype_id").iloc[:, 3:]
            cov[r] = pd.read_csv(d / "covariates.tsv", sep="\t", index_col=0).T
        for a, b in PAIRS:
            if args.at:
                e = pd.read_csv(ROOT / f"results/06-compare/eqtl-comparison/{args.at}/paired-{args.variants}/"
                                f"robust-{a}-vs-{b}.tsv", sep="\t", dtype={"gene_ncbi": str})
                e = e[~e.artefact_flag.astype(bool)].rename(columns={"z_paired": "disc_z", "delta": "disc_delta",
                                                                      "robust_indep": "in_indep_set"})
            else:
                e = pd.read_csv(base / f"{m}-{args.variants}/effects-{a}-vs-{b}.tsv", sep="\t",
                                dtype={"gene_ncbi": str})
            e["pheno_a"], e["pheno_b"] = e.gene_ncbi.map(gid[a]).values, e.gene_ncbi.map(gid[b]).values
            e = e[e.pheno_a.isin(ph[a].index) & e.pheno_b.isin(ph[b].index)
                  & e.var_a.isin(rd[GENO[a]].variant_df.index) & e.var_b.isin(rd[GENO[b]].variant_df.index)].copy()
            samples = [s for s in ph[a].columns if s in ph[b].columns and s in cov[a].index and s in cov[b].index
                       and s in rd["t2t"].sample_ids and s in rd["grch38"].sample_ids and s not in shared]
            if cnt is not None:
                e = e[e.symbol.isin(cnt.symbol)].reset_index(drop=True)
                cn_of = {s: x[samples].to_numpy(dtype=float).T for s, x in cnt.groupby("symbol")}
            Q = {}
            for r in (a, b):
                Q[r], _ = np.linalg.qr(np.column_stack([np.ones(len(samples)), cov[r].loc[samples].values]))
            dof = len(samples) - Q[a].shape[1] - 1  # covariate sets have the same size on every reference
            G = {}
            for r, v in ((a, e.var_a), (b, e.var_b)):
                s = GENO[r]
                u = v.unique().tolist()
                x = pd.DataFrame(rd[s].read_list(u, impute_mean=True), index=u, columns=rd[s].sample_ids)
                G[r] = resid(x.loc[v.values, samples].values.astype(float), Q[r])
            G[b] = G[b] * np.sign(e.geno_r.values)[:, None]
            Y = {a: resid(ph[a].loc[e.pheno_a.values, samples].values.astype(float), Q[a]),
                 b: resid(ph[b].loc[e.pheno_b.values, samples].values.astype(float), Q[b])}
            dof = np.full(len(e), dof)
            if cnt is not None:
                cn_r2 = np.empty(len(e))
                for i, s in enumerate(e.symbol):
                    for r in (a, b):
                        K, _ = np.linalg.qr(resid(cn_of[s].T, Q[r]).T)
                        K = K[:, np.abs(np.diag(_)) > 1e-8]  # drop a family column collinear with the gene column
                        if r == a:
                            g0 = (G[a][i] ** 2).sum()
                        G[r][i] -= K @ (K.T @ G[r][i])
                        Y[r][i] -= K @ (K.T @ Y[r][i])
                    cn_r2[i] = 1 - (G[a][i] ** 2).sum() / g0
                    dof[i] -= K.shape[1]
            S = {r: (G[r] ** 2).sum(1) for r in (a, b)}
            sl = {r: (G[r] * Y[r]).sum(1) / S[r] for r in (a, b)}
            E = {r: Y[r] - sl[r][:, None] * G[r] for r in (a, b)}
            se = {r: np.sqrt((E[r] ** 2).sum(1) / dof / S[r]) for r in (a, b)}
            c = (E[a] * E[b]).sum(1) / dof * (G[a] * G[b]).sum(1) / (S[a] * S[b])
            x = e[["gene_ncbi", "symbol", "biotype", "region_class", "artefact_flag", "var_a", "var_b", "geno_r",
                   "slope_a", "slope_b", "z_diff", "p_diff", "changed"]
                  + (["disc_z", "disc_delta", "in_indep_set"] if args.at else [])].copy()
            chk = max(np.abs(sl[a] - x.slope_a.values).max(), np.abs(sl[b] - x.slope_b.values).max())
            x["rho"] = (E[a] * E[b]).sum(1) / np.sqrt((E[a] ** 2).sum(1) * (E[b] ** 2).sum(1))
            if cnt is not None:
                x["slope_a"], x["slope_b"], x["cn_r2"] = sl[a], sl[b], cn_r2
            if args.at:  # discovery statistics replaced by this cohort's
                x["slope_a"], x["slope_b"] = sl[a], sl[b]
                x["z_diff"] = (sl[a] - sl[b]) / np.sqrt(se[a] ** 2 + se[b] ** 2)
                x["p_diff"] = 2 * norm.sf(x.z_diff.abs())
                x["changed"] = x.p_diff * len(x) < 0.05
            x["se_a"], x["se_b"], x["cov_ab"] = se[a], se[b], c
            x["se_diff_indep"] = np.sqrt(se[a] ** 2 + se[b] ** 2)
            x["se_diff_paired"] = np.sqrt(np.maximum(se[a] ** 2 + se[b] ** 2 - 2 * c, 1e-12))
            x["delta"] = x.slope_a - x.slope_b
            x["z_paired"] = x.delta / x.se_diff_paired
            x["p_paired"] = 2 * norm.sf(x.z_paired.abs())
            x["changed_paired"] = x.p_paired * len(x) < 0.05
            x.sort_values("p_paired").to_csv(out / f"effects-{a}-vs-{b}-{m}.tsv", sep="\t", index=False,
                                             float_format="%.4g")
            res[(m, a, b)] = x.set_index("gene_ncbi")
            lam = lambda z: np.median(z ** 2) / 0.4549
            print(f"{m} {a} vs {b}: {len(x)} genes, {len(samples)} individuals; slope check max diff {chk:.2g}; "
                  f"median rho {np.median(x.rho):.3f}; median se ratio paired/indep "
                  f"{np.median(x.se_diff_paired / x.se_diff_indep):.3f}; lambda indep {lam(x.z_diff):.2f}, "
                  f"paired {lam(x.z_paired):.2f}; changed indep {int(x.changed.sum())}, paired "
                  f"{int(x.changed_paired.sum())}", flush=True)

            # bootstrap over individuals for a random subset of genes
            k = rng.choice(len(x), size=min(args.nboot_genes, len(x)), replace=False)
            D = np.empty((len(k), args.nboot))
            for i in range(args.nboot):
                j = rng.integers(0, len(samples), len(samples))
                ga, gb, ya, yb = G[a][k][:, j], G[b][k][:, j], Y[a][k][:, j], Y[b][k][:, j]
                D[:, i] = (ga * ya).sum(1) / (ga ** 2).sum(1) - (gb * yb).sum(1) / (gb ** 2).sum(1)
            sd = D.std(1, ddof=1)
            calib.append({"method": m, "pair": f"{a} vs {b}", "genes": len(k), "boot": args.nboot,
                          "median_boot_over_paired": np.median(sd / x.se_diff_paired.values[k]),
                          "median_boot_over_indep": np.median(sd / x.se_diff_indep.values[k]),
                          "q10_boot_over_paired": np.quantile(sd / x.se_diff_paired.values[k], 0.1),
                          "q90_boot_over_paired": np.quantile(sd / x.se_diff_paired.values[k], 0.9),
                          "lambda_indep": lam(x.z_diff), "lambda_paired": lam(x.z_paired),
                          "median_rho": np.median(x.rho)})
    pd.DataFrame(calib).to_csv(out / "calibration.tsv", sep="\t", index=False, float_format="%.4g")

    if args.at:
        replication(res, out)
        return
    rows = []
    for a, b in PAIRS:
        s, q = res[("star", a, b)], res[("salmon", a, b)]
        j = s.join(q[["delta", "z_paired", "changed_paired", "changed"]], rsuffix="_salmon", how="inner")
        same = np.sign(j.z_paired) == np.sign(j.z_paired_salmon)
        j["robust_paired"] = j.changed_paired & j.changed_paired_salmon & same
        j["robust_indep"] = j.changed & j.changed_salmon & (np.sign(j.z_diff) == np.sign(j.z_paired_salmon))
        r = j[j.robust_paired].sort_values("p_paired")
        r.to_csv(out / f"robust-{a}-vs-{b}.tsv", sep="\t", float_format="%.4g")
        u = r[~r.artefact_flag.astype(bool)]
        ad = u.delta.abs()
        rows.append({"pair": f"{a} vs {b}", "genes_both_methods": len(j), "robust_indep": int(j.robust_indep.sum()),
                     "robust_paired": len(r), "robust_paired_unflagged": len(u),
                     "unflagged_also_indep": int(u.robust_indep.sum()),
                     "unflagged_new": int((~u.robust_indep).sum()),
                     "unflagged_abs_delta_ge_0.1": int((ad >= 0.1).sum()),
                     "unflagged_abs_delta_ge_0.2": int((ad >= 0.2).sum()),
                     "median_abs_delta_new": float(u.loc[~u.robust_indep, "delta"].abs().median())})
    s = pd.DataFrame(rows)
    s.to_csv(out / "summary.tsv", sep="\t", index=False, float_format="%.4g")
    pd.set_option("display.width", 250)
    print("\n" + pd.DataFrame(calib).round(3).to_string(index=False))
    print("\n" + s.to_string(index=False))


if __name__ == "__main__":
    main()
