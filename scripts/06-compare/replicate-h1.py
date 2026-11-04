"""Step 6: replication of the H1 genes (lead variant in T2T-only sequence) of a discovery cohort, at the same variants.

    pixi run python scripts/06-compare/replicate-h1.py [--discovery mage] [--replication geuvadis] [--method star]

Discovery: results/06-compare/eqtl-comparison/<discovery>/h1-variant-qc.tsv (h1-variant-qc.py): per gene the T2T-only
lead variant, its slope, the best shared (matched) variant, the genotype call and the signal call.
In the replication cohort (T2T, all variants), for exactly the same variants (same IDs and alleles, one 1KGP panel):
  rep_slope, rep_p          nominal association of the discovery lead variant (map_nominal parquet)
  rep_shared_p              nominal p of the discovery best shared variant
  rep_lead_p_given_shared   lead p after conditioning on that shared variant (covariate-residualised OLS, as in
                            h1-variant-qc.py)
  rep_het_F, rep_hwe_p      heterozygosity and within-population HWE of the lead in the replication individuals
  rep_qval, rep_lead        the gene's permutation q and its own lead variant in the replication cohort
Calls, with n = genes tested:
  replicated        same sign as discovery and rep_p < 0.05 / n
  signal_replicated replicated and rep_lead_p_given_shared < 0.05 / n (or no shared variant / perfect proxy)
Writes results/06-compare/eqtl-comparison/replication-<replication>/h1-genes.tsv and prints a summary.
--exclude-shared drops the discovery cohort's individuals and recomputes rep_slope, rep_p and rep_shared_p on the rest
(replicate-changes.py ols_nominal); output in replication-<replication>-independent/.
"""
import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from tensorqtl import pgen

ROOT = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location("h1qc", ROOT / "scripts/06-compare/h1-variant-qc.py")
h1qc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(h1qc)
_spec = importlib.util.spec_from_file_location("repch", ROOT / "scripts/06-compare/replicate-changes.py")
repch = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(repch)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--discovery", default="mage")
    ap.add_argument("--replication", default="geuvadis")
    ap.add_argument("--method", default="star")
    ap.add_argument("--exclude-shared", action="store_true", help="drop the discovery cohort's individuals")
    args = ap.parse_args()
    base = ROOT / "results/06-compare/eqtl-comparison"
    out = base / f"replication-{args.replication}{'-independent' if args.exclude_shared else ''}"
    exclude = repch.shared_individuals(args.discovery, args.replication) if args.exclude_shared else set()
    out.mkdir(parents=True, exist_ok=True)
    disc = pd.read_csv(base / f"{args.discovery}/h1-variant-qc.tsv", sep="\t")

    d = ROOT / f"results/05-eqtl/{args.replication}/chm13v2.0/{args.method}"
    runs = pd.read_csv(ROOT / f"results/04-expression/{args.replication}/matrices/runs.tsv", sep="\t")
    pop = runs.drop_duplicates("sample").set_index("sample")["ancestry"].fillna("unknown").str.split(" ").str[0]
    ph = pd.read_csv(d / "phenotypes.bed.gz", sep="\t").set_index("phenotype_id")
    cov = pd.read_csv(d / "covariates.tsv", sep="\t", index_col=0).T
    rd = pgen.PgenReader(str(ROOT / f"results/05-eqtl/{args.replication}/genotypes/t2t.all"))
    samples = [s for s in rd.sample_ids if s in cov.index and s in ph.columns and s not in exclude]
    print(f"{len(samples)} {args.replication} individuals ({len(exclude)} shared with {args.discovery} excluded)")
    sidx = [rd.sample_ids.index(s) for s in samples]
    C = np.column_stack([np.ones(len(samples)), cov.loc[samples].values])
    pops = pop.reindex(samples)
    perm = pd.read_csv(d / "all/cis-eqtl.tsv.gz", sep="\t").set_index("phenotype_id")
    have = set(rd.variant_ids)

    def geno(v):
        g = np.asarray(rd.read(v), dtype=float)[sidx]
        return g

    nominal = {}
    rows = []
    for r in disc.itertuples():
        chrom = r.variant_id.split(":")[0]
        rec = {"gene": r.gene, "variant_id": r.variant_id, "disc_slope": r.slope, "disc_genotype_call": r.genotype_call,
               "disc_signal": r.signal, "best_shared_variant": r.best_shared_variant}
        if r.gene in perm.index:
            rec["rep_qval"], rec["rep_lead"] = perm.loc[r.gene, "qval"], perm.loc[r.gene, "variant_id"]
        if r.gene not in ph.index or r.variant_id not in have:
            rows.append(rec)  # gene not tested or variant below MAF 0.01 in the replication cohort
            continue
        if chrom not in nominal:
            nominal[chrom] = pd.read_parquet(d / f"all/nominal/{chrom}.cis_qtl_pairs.{chrom}.parquet",
                                             columns=["phenotype_id", "variant_id", "pval_nominal", "slope",
                                                      "slope_se"],
                                             filters=[("phenotype_id", "in", list(disc.gene))]
                                             ).set_index(["phenotype_id", "variant_id"])
        n = nominal[chrom]
        if args.exclude_shared:
            w = pd.DataFrame({"phenotype_id": r.gene, "variant_id": [v for v in (r.variant_id, r.best_shared_variant)
                                                                      if isinstance(v, str)]})
            n = repch.ols_nominal(args.replication, "chm13v2.0", args.method, "all", w, exclude)
        if (r.gene, r.variant_id) not in n.index:
            rows.append(rec)
            continue
        x = n.loc[(r.gene, r.variant_id)]
        rec.update(rep_slope=x.slope, rep_slope_se=x.slope_se, rep_p=x.pval_nominal)
        g = geno(r.variant_id)
        rec["rep_het_F"], rec["rep_hwe_p"] = h1qc.het_test(g, pops)
        s = r.best_shared_variant
        if isinstance(s, str) and (r.gene, s) in n.index:
            rec["rep_shared_p"] = n.loc[(r.gene, s), "pval_nominal"]
            y = h1qc.residualise(ph.loc[r.gene, samples].values.astype(float), C)
            gl = h1qc.residualise(np.nan_to_num(g, nan=np.nanmean(g)), C)
            gs = geno(s)
            gs = h1qc.residualise(np.nan_to_num(gs, nan=np.nanmean(gs)), C)
            rec["rep_ld_r2"] = np.corrcoef(gl, gs)[0, 1] ** 2
            if rec["rep_ld_r2"] <= 0.999:
                X = np.column_stack([gl, gs])
                beta, *_ = np.linalg.lstsq(X, y, rcond=None)
                res = y - X @ beta
                dof = len(y) - C.shape[1] - 2
                se = np.sqrt(res @ res / dof * np.linalg.inv(X.T @ X).diagonal())
                rec["rep_lead_p_given_shared"] = 2 * stats.t.sf(abs(beta[0] / se[0]), dof)
                rec["rep_shared_p_given_lead"] = 2 * stats.t.sf(abs(beta[1] / se[1]), dof)
        rows.append(rec)

    t = pd.DataFrame(rows)
    for c in ("rep_p", "rep_slope", "rep_lead_p_given_shared", "rep_ld_r2", "rep_shared_p"):
        if c not in t:
            t[c] = np.nan
    t["tested"] = t.rep_p.notna()
    k = int(t.tested.sum())
    t["same_sign"] = t.tested & (np.sign(t.rep_slope) == np.sign(t.disc_slope))
    t["replicated"] = t.same_sign & (t.rep_p < 0.05 / max(k, 1))
    no_sep = t.best_shared_variant.isna() | t.rep_shared_p.isna() | (t.rep_ld_r2 > 0.999)
    t["signal_replicated"] = t.replicated & (no_sep | (t.rep_lead_p_given_shared < 0.05 / max(k, 1)))
    t["rep_lead_is_same"] = t.rep_lead.eq(t.variant_id)
    t.to_csv(out / "h1-genes.tsv", sep="\t", index=False, float_format="%.4g")

    pd.set_option("display.width", 250)
    print(f"{len(t)} discovery H1 genes, {k} tested in {args.replication} (Bonferroni 0.05/{k})")
    agg = {"genes": ("gene", "size"), "tested": ("tested", "sum"), "same_sign": ("same_sign", "sum"),
           "replicated": ("replicated", "sum"), "signal_replicated": ("signal_replicated", "sum"),
           "same_lead": ("rep_lead_is_same", "sum")}
    print(t.groupby("disc_genotype_call").agg(**agg).to_string())
    solid = t[(t.disc_genotype_call == "solid") & t.disc_signal.isin(["T2T-only variant explains it", "both contribute"])]
    print(f"\nSolid genes where the T2T-only variant carries the signal in {args.discovery}: {len(solid)}")
    print(solid.groupby(lambda i: "all").agg(**agg).to_string(index=False))
    print(solid[["gene", "variant_id", "disc_slope", "rep_slope", "rep_p", "rep_lead_p_given_shared", "rep_het_F",
                 "rep_qval", "rep_lead_is_same", "replicated", "signal_replicated"]].to_string(index=False))


if __name__ == "__main__":
    main()
