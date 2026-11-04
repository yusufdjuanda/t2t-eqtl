"""Step 6 (H6): colocalisation of the changed genes' eQTLs with GWAS, on each reference.

    pixi run python scripts/06-compare/coloc.py [--cohort mage] [--pair chm13v2.0 grch38-ffx]

Question: does the gene a GWAS locus points to, or the direction of that link, change with the reference?
Method: coloc's approximate Bayes factor test (coloc.abf; Giambartolomei et al. 2014, Wakefield 2009), one causal
variant per trait in the region, written out here:
  per variant and trait, lABF = 0.5 log(1 - r) + 0.5 r z^2 with r = W / (W + se^2);
  prior effect SD: 0.15 for the eQTL (expression is inverse-normal, SD 1), 0.2 for GWAS log odds ratios (coloc defaults);
  priors p1 = p2 = 1e-4, p12 = 1e-5 (coloc defaults); posteriors PP0-PP4, PP4 = one variant shared by both.
Variants: the gene's whole cis window (tensorQTL map_nominal, STAR, matched variants), as GRCh38 IDs; T2T IDs are
converted with the variant matching table, so both references are tested on the same variants. GWAS variants are
joined on position and alleles (either orientation).
Direction: at the variant with the largest shared evidence (lABF_eQTL + lABF_GWAS), sign(eQTL slope x GWAS beta) with
both on the same allele: + means the allele that raises expression raises risk. Allele orientation of the T2T slope is
taken from the matching table (allele_match) and reported for SNPs only.
Genes: all tiers of gene-sets/changed-<A>-vs-<B>.tsv. GWAS: data/external/gwas-sumstats/<GCST>/*.h.tsv.gz (TRAITS).
Writes results/06-compare/gwas/coloc-<cohort>-<A>-vs-<B>.tsv (gene x trait, both references side by side).
"""
import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[2]
TRAITS = {"GCST90627749": "hypothyroidism", "GCST003156": "systemic lupus erythematosus",
          "GCST004133": "ulcerative colitis", "GCST90558092": "multiple sclerosis",
          "GCST006867": "type 2 diabetes", "GCST90624412": "atrial fibrillation"}
W_EQTL, W_GWAS = 0.15 ** 2, 0.2 ** 2
P1, P2, P12 = 1e-4, 1e-4, 1e-5


def labf(z, se, w):
    r = w / (w + se ** 2)
    return 0.5 * np.log(1 - r) + 0.5 * r * z ** 2


def lsum(x):
    m = np.max(x)
    return m + np.log(np.sum(np.exp(x - m)))


def coloc_abf(l1, l2):
    """Posterior probabilities PP0-PP4 from the per-variant log Bayes factors of two traits."""
    s1, s2, s12 = lsum(l1), lsum(l2), lsum(l1 + l2)
    h = np.array([0.0, np.log(P1) + s1, np.log(P2) + s2,
                  np.log(P1) + np.log(P2) + s1 + s2 + np.log1p(-min(np.exp(s12 - s1 - s2), 1 - 1e-16)),
                  np.log(P12) + s12])
    return np.exp(h - lsum(h))


def read_gwas(path, windows):
    """Harmonised GWAS-SSF file -> chrom, pos, a_other, a_effect, beta, se for variants inside the windows."""
    cols = pd.read_csv(path, sep="\t", nrows=0).columns
    pick = lambda *c: next(x for x in c if x in cols)
    cmap = {pick("hm_chrom", "chromosome"): "chrom", pick("hm_pos", "base_pair_location"): "pos",
            pick("hm_other_allele", "other_allele"): "a_other", pick("hm_effect_allele", "effect_allele"): "a_effect",
            pick("hm_beta", "beta"): "beta", pick("standard_error"): "se", pick("p_value"): "p"}
    if "hm_odds_ratio" in cols or "odds_ratio" in cols:
        cmap[pick("hm_odds_ratio", "odds_ratio")] = "or"
    parts = []
    for ch in pd.read_csv(path, sep="\t", usecols=list(cmap), dtype=str, chunksize=2_000_000):
        ch = ch.rename(columns=cmap)
        ch["chrom"] = "chr" + ch.chrom.str.replace("chr", "", regex=False)
        ch["pos"] = pd.to_numeric(ch.pos, errors="coerce")
        keep = np.zeros(len(ch), bool)
        for c, iv in windows.items():
            s = (ch.chrom == c).values
            if not s.any():
                continue
            p = ch.pos.values[s]
            k = np.zeros(s.sum(), bool)
            for lo, hi in iv:
                k |= (p >= lo) & (p <= hi)
            keep[np.flatnonzero(s)[k]] = True
        parts.append(ch[keep])
    g = pd.concat(parts)
    for c in ("beta", "se", "p", "or"):
        if c in g:
            g[c] = pd.to_numeric(g[c], errors="coerce")
    if g.beta.isna().all() and "or" in g:
        g["beta"] = np.log(g["or"])
    miss = g.se.isna() | (g.se <= 0)
    if miss.any():  # se from beta and p
        z = norm.isf(g.loc[miss, "p"] / 2)
        g.loc[miss, "se"] = (g.loc[miss, "beta"].abs() / z).replace(0, np.nan)
    g = g.dropna(subset=["pos", "beta", "se"])
    g = g[g.a_other.str.fullmatch("[ACGT]+") & g.a_effect.str.fullmatch("[ACGT]+")]
    g["pos"] = g.pos.astype(int)
    return g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", default="mage")
    ap.add_argument("--pair", nargs=2, default=["chm13v2.0", "grch38-ffx"])
    ap.add_argument("--traits", nargs="+", default=list(TRAITS), help="GCST accessions (default: all of TRAITS)")
    ap.add_argument("--out-suffix", default="")
    args = ap.parse_args()
    a, b = args.pair
    out = ROOT / "results/06-compare/gwas"
    out.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(ROOT / f"results/06-compare/eqtl-comparison/{args.cohort}/gene-sets/changed-{a}-vs-{b}.tsv",
                        sep="\t", dtype={"gene_ncbi": str})

    m = pd.concat(pd.read_csv(f, sep="\t", usecols=["side", "category", "t2t_id", "grch38_id", "allele_match"],
                              low_memory=False)
                  for f in sorted((ROOT / "results/06-compare/variant-matching").glob("chr*.tsv.gz")))
    m = m[(m.side == "t2t") & (m.category == "matched")]
    t2g = dict(zip(m.t2t_id, m.grch38_id))
    t2flip = dict(zip(m.t2t_id, np.where(m.allele_match.eq("swapped"), -1.0, 1.0)))
    del m

    # eQTL statistics per reference, as GRCh38 IDs
    eq = {}
    for ref in (a, b):
        gt = pd.read_csv(ROOT / f"data/star-references/{ref}/genes.tsv", sep="\t", dtype={"gene_ncbi": str})
        pid = gt[gt.gene_ncbi.isin(genes.gene_ncbi)].set_index("gene_id").gene_ncbi
        parts = []
        for f in sorted(glob.glob(str(ROOT / f"results/05-eqtl/{args.cohort}/{ref}/star/matched/nominal/*.parquet"))):
            x = pd.read_parquet(f, columns=["phenotype_id", "variant_id", "slope", "slope_se"],
                                filters=[("phenotype_id", "in", list(pid.index))])
            parts.append(x)
        x = pd.concat(parts)
        x["gene_ncbi"] = x.phenotype_id.map(pid)
        if ref == "chm13v2.0":
            x["flip"] = x.variant_id.map(t2flip)
            x["variant_id"] = x.variant_id.map(t2g)
            x["slope"] = x.slope * x.flip
            x = x.dropna(subset=["variant_id"])
        eq[ref] = x[["gene_ncbi", "variant_id", "slope", "slope_se"]]
        print(f"{ref}: {x.gene_ncbi.nunique()} genes, {len(x)} gene-variant pairs", flush=True)

    allv = pd.concat(eq.values())
    p = allv.variant_id.str.split(":", expand=True)
    allv = allv.assign(chrom=p[0], pos=p[1].astype(int))
    windows = {c: x.groupby("gene_ncbi").pos.agg(["min", "max"]).values.tolist() for c, x in allv.groupby("chrom")}

    rows = []
    for gcst in args.traits:
        trait = TRAITS[gcst]
        f = glob.glob(str(ROOT / f"data/external/gwas-sumstats/{gcst}/*.h.tsv.gz"))
        if not f:
            print(f"{gcst} ({trait}): no file, skipped", flush=True)
            continue
        gw = read_gwas(f[0], windows)
        fwd = gw.chrom + ":" + gw.pos.astype(str) + ":" + gw.a_other + ":" + gw.a_effect
        rev = gw.chrom + ":" + gw.pos.astype(str) + ":" + gw.a_effect + ":" + gw.a_other
        gb = pd.concat([pd.DataFrame({"variant_id": fwd, "g_beta": gw.beta, "g_se": gw.se}),
                        pd.DataFrame({"variant_id": rev, "g_beta": -gw.beta, "g_se": gw.se})]
                       ).drop_duplicates("variant_id").set_index("variant_id")
        print(f"{gcst} ({trait}): {len(gw)} GWAS variants in the windows", flush=True)
        for gene, _ in genes.groupby("gene_ncbi"):
            rec = {"gene_ncbi": gene, "gcst": gcst, "trait": trait}
            for ref in (a, b):
                e = eq[ref][eq[ref].gene_ncbi == gene].set_index("variant_id")
                j = e.join(gb, how="inner")
                rec[f"nsnps_{ref}"] = len(j)
                if len(j) < 50:
                    continue
                l1 = labf(j.slope / j.slope_se, j.slope_se, W_EQTL).values
                l2 = labf(j.g_beta / j.g_se, j.g_se, W_GWAS).values
                pp = coloc_abf(l1, l2)
                for k in range(5):
                    rec[f"PP{k}_{ref}"] = pp[k]
                i = int(np.argmax(l1 + l2))
                v = j.index[i]
                snp = all(len(x) == 1 for x in v.split(":")[2:])
                rec[f"top_{ref}"] = v
                rec[f"eqtl_z_{ref}"] = j.slope.iloc[i] / j.slope_se.iloc[i]
                rec[f"gwas_z_{ref}"] = j.g_beta.iloc[i] / j.g_se.iloc[i]
                rec[f"direction_{ref}"] = (np.sign(j.slope.iloc[i] * j.g_beta.iloc[i]) if snp or ref != "chm13v2.0"
                                           else np.nan)
                rec[f"gwas_min_p_{ref}"] = 2 * norm.sf(np.abs(j.g_beta / j.g_se)).min()
            rows.append(rec)

    t = pd.DataFrame(rows).merge(genes[["gene_ncbi", "symbol", "tier", "region_class", "slope_a", "slope_b"]],
                                 on="gene_ncbi")
    for ref in (a, b):
        if f"PP4_{ref}" not in t:
            t[f"PP4_{ref}"] = np.nan
    pa, pb = t[f"PP4_{a}"].fillna(0), t[f"PP4_{b}"].fillna(0)
    t["call"] = np.select([(pa >= 0.8) & (pb >= 0.8), (pa >= 0.8) & (pb < 0.5), (pb >= 0.8) & (pa < 0.5),
                           (pa >= 0.5) | (pb >= 0.5)],
                          ["both", f"{a} only", f"{b} only", "suggestive"], "none")
    t["direction_flips"] = (t.call == "both") & (t[f"direction_{a}"] * t[f"direction_{b}"] < 0)
    t.to_csv(out / f"coloc-{args.cohort}-{a}-vs-{b}{args.out_suffix}.tsv", sep="\t", index=False, float_format="%.4g")

    pd.set_option("display.width", 250)
    print("\ncalls (PP4 >= 0.8; 'only' = < 0.5 on the other reference), genes x traits:")
    print(pd.crosstab(t.trait, t.call).to_string())
    for tiers, lab in ((["core", "main"], "main"), (["core", "main", "extent"], "extent")):
        x = t[t.tier.isin(tiers)]
        print(f"{lab}: genes colocalising (PP4 >= 0.8) on either reference: "
              f"{x.loc[x.call.isin(['both', f'{a} only', f'{b} only']), 'gene_ncbi'].nunique()}; "
              f"on one reference only: {x.loc[x.call.isin([f'{a} only', f'{b} only']), 'gene_ncbi'].nunique()}; "
              f"direction flips: {x.loc[x.direction_flips, 'gene_ncbi'].nunique()}")
    x = t[t.call != "none"].sort_values(["call", "trait"])
    print(x[["symbol", "tier", "region_class", "trait", "call", f"PP4_{a}", f"PP4_{b}", f"direction_{a}",
             f"direction_{b}", f"nsnps_{a}", f"gwas_min_p_{a}"]].to_string(index=False))


if __name__ == "__main__":
    main()
