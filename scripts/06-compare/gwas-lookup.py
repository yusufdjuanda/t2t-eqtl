"""Step 6: GWAS lookup for the changed genes: are their eQTL variants, or variants in LD with them, GWAS hits?

    pixi run python scripts/06-compare/gwas-lookup.py [--cohort mage] [--pair chm13v2.0 grch38-ffx]

A first screen before colocalisation (H6): it shows which changed genes sit at known disease or trait loci, and for
which traits, so that full summary statistics can be fetched for those. Being in LD with a GWAS hit does not show a
shared causal variant; coloc does.
Genes: all tiers of gene-sets/changed-<A>-vs-<B>.tsv (define-gene-sets.py).
Query variants per gene, all as GRCh38 IDs (the two GRCh38 references share coordinates; T2T IDs are converted with
the variant matching table, matched variants only):
  pair_a, pair_b    the variant pair compared (paired-test.py)
  lead_<ref>        each reference's own lead variant (STAR, matched variants, map_cis)
GWAS: GWAS Catalog associations (data/external/gwas-catalog/, release 'latest' fetched 2026-10-29), GRCh38
positions, single-SNP rows. For every association within 1 Mb of a query variant, r^2 between the two in the cohort's
GRCh38 genotypes (all individuals, and the EUR subset, as most GWAS are European); a proxy hit is r^2 >= 0.8 in either.
--background K: the same for K unchanged genes per core/main gene (same region class and eQTL-strength quintile),
the comparison rate. Writes results/06-compare/gwas/lookup-<cohort>-<A>-vs-<B>[-background].tsv (gene x association, r^2 >= 0.5) and
genes-<cohort>-<A>-vs-<B>.tsv (one row per gene), and prints a summary.
"""
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pgenlib

ROOT = Path(__file__).resolve().parents[2]
WINDOW = 1_000_000
IMMUNE = re.compile(r"lupus|rheumatoid|multiple sclerosis|crohn|ulcerative colitis|inflammatory bowel|type 1 diabetes|"
                    r"asthma|psoriasis|celiac|coeliac|ankylosing|vitiligo|graves|hashimoto|hypothyroid|primary biliary|"
                    r"sclerosing cholangitis|allerg|eczema|atopic|sjogren|myasthenia|autoimmun|lymphocyte|leukocyte|"
                    r"white blood cell|neutrophil|eosinophil|basophil|monocyte|immunoglobulin|lymphoma|leukemia|"
                    r"leukaemia|kawasaki|sarcoidosis|spondyl|juvenile idiopathic", re.I)


def background(args, a, b, changed, seed=1):
    """K unchanged compared genes per core/main gene, drawn from the same region class and quintile of eQTL strength
    (max -log10 p of the compared pair), from the STAR effects table; artefact-flagged and extent genes excluded."""
    e = pd.read_csv(ROOT / f"results/06-compare/eqtl-comparison/{args.cohort}/star-matched/effects-{a}-vs-{b}.tsv",
                    sep="\t", dtype={"gene_ncbi": str})
    e = e[~e.artefact_flag.astype(bool) & e.var_a.notna() & e.var_b.notna()].copy()
    e["strength"] = e[["log10p_a", "log10p_b"]].max(axis=1)
    edges = np.quantile(e.strength, [0.2, 0.4, 0.6, 0.8])
    e["bin"] = np.digitize(e.strength, edges)
    target = e[e.gene_ncbi.isin(changed.loc[changed.tier.isin(["core", "main"]), "gene_ncbi"])]
    pool = e[~e.gene_ncbi.isin(changed.gene_ncbi)]
    rng = np.random.default_rng(seed)
    pick = []
    for (rc, bn), n in target.groupby(["region_class", "bin"]).size().items():
        cand = pool[(pool.region_class == rc) & (pool.bin == bn)]
        pick.append(cand.iloc[rng.choice(len(cand), size=min(len(cand), n * args.background), replace=False)])
    bg = pd.concat(pick).assign(tier="background")
    print(f"background: {len(bg)} unchanged genes for {len(target)} core/main genes; region classes "
          f"{bg.region_class.value_counts().to_dict()}", flush=True)
    return bg[["gene_ncbi", "symbol", "biotype", "region_class", "tier", "var_a", "var_b", "slope_a", "slope_b"]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", default="mage")
    ap.add_argument("--pair", nargs=2, default=["chm13v2.0", "grch38-ffx"])
    ap.add_argument("--background", type=int, default=0,
                    help="instead of the changed genes, K unchanged genes per core/main gene (matched on region class "
                         "and eQTL strength); output suffix -background")
    args = ap.parse_args()
    a, b = args.pair
    out = ROOT / "results/06-compare/gwas"
    out.mkdir(parents=True, exist_ok=True)

    g = pd.read_csv(ROOT / f"results/06-compare/eqtl-comparison/{args.cohort}/gene-sets/changed-{a}-vs-{b}.tsv",
                    sep="\t", dtype={"gene_ncbi": str})
    tag = ""
    if args.background:
        g = background(args, a, b, g)
        tag = "-background"
    m = pd.concat(pd.read_csv(f, sep="\t", usecols=["side", "category", "t2t_id", "grch38_id"], low_memory=False)
                  for f in sorted((ROOT / "results/06-compare/variant-matching").glob("chr*.tsv.gz")))
    t2g = m[(m.side == "t2t") & (m.category == "matched")].set_index("t2t_id").grch38_id.to_dict()
    to38 = {r: (lambda v: t2g.get(v)) if r == "chm13v2.0" else (lambda v: v) for r in (a, b)}

    q = []
    for r in g.itertuples():
        q.append((r.gene_ncbi, "pair_a", to38[a](r.var_a)))
        q.append((r.gene_ncbi, "pair_b", to38[b](r.var_b)))
    for ref in (a, b):
        genes = pd.read_csv(ROOT / f"data/star-references/{ref}/genes.tsv", sep="\t", dtype={"gene_ncbi": str})
        gid = genes.set_index("gene_id").gene_ncbi
        lead = pd.read_csv(ROOT / f"results/05-eqtl/{args.cohort}/{ref}/star/matched/cis-eqtl.tsv.gz", sep="\t")
        lead["gene_ncbi"] = lead.phenotype_id.map(gid)
        lead = lead[lead.gene_ncbi.isin(g.gene_ncbi)]
        q += [(x.gene_ncbi, f"lead_{ref}", to38[ref](x.variant_id)) for x in lead.itertuples()]
    q = pd.DataFrame(q, columns=["gene_ncbi", "role", "query"]).dropna().drop_duplicates()
    print(f"{len(g)} genes, {len(q)} query variants ({q['query'].nunique()} distinct)", flush=True)

    # GWAS Catalog, single-SNP rows with a GRCh38 position
    gw = pd.read_csv(ROOT / "data/external/gwas-catalog/gwas-catalog-download-associations-alt-full.tsv", sep="\t",
                     usecols=["PUBMEDID", "FIRST AUTHOR", "DISEASE/TRAIT", "MAPPED_TRAIT", "CHR_ID", "CHR_POS", "SNPS",
                              "STRONGEST SNP-RISK ALLELE", "P-VALUE", "PVALUE_MLOG", "OR or BETA", "STUDY ACCESSION"],
                     dtype=str, low_memory=False).rename(columns={
        "PUBMEDID": "pubmed", "FIRST AUTHOR": "author", "DISEASE/TRAIT": "trait", "MAPPED_TRAIT": "mapped_trait",
        "SNPS": "snp", "STRONGEST SNP-RISK ALLELE": "risk_allele", "OR or BETA": "or_beta", "STUDY ACCESSION": "study"})
    gw = gw[gw.CHR_ID.str.fullmatch(r"\d+|X", na=False) & gw.CHR_POS.str.fullmatch(r"\d+", na=False)].copy()
    gw["chrom"], gw["pos"] = "chr" + gw.CHR_ID, gw.CHR_POS.astype(int)
    gw["mlog10p"] = pd.to_numeric(gw.PVALUE_MLOG, errors="coerce")
    print(f"GWAS Catalog: {len(gw)} single-SNP associations with a GRCh38 position", flush=True)

    # cohort GRCh38 genotypes
    gdir = ROOT / f"results/05-eqtl/{args.cohort}/genotypes"
    pv = pd.read_csv(gdir / "grch38.all.pvar", sep="\t", comment="#", header=None, usecols=[0, 1, 2],
                     names=["chrom", "pos", "id"], dtype={"chrom": str})
    pv["idx"] = np.arange(len(pv))
    id2idx = dict(zip(pv.id, pv.idx))
    psam = pd.read_csv(gdir / "grch38.all.psam", sep="\t")
    runs = pd.read_csv(ROOT / f"results/04-expression/{args.cohort}/matrices/runs.tsv", sep="\t")
    eur_ids = set(runs.loc[runs.ancestry.fillna("").str.contains("EUR"), "sample"])
    eur = psam.iloc[:, 0].isin(eur_ids).values
    print(f"{len(psam)} individuals, {int(eur.sum())} EUR", flush=True)
    rd = pgenlib.PgenReader(str(gdir / "grch38.all.pgen").encode())
    n = rd.get_raw_sample_ct()
    cache = {}

    def dos(i):
        if i not in cache:
            d = np.empty(n, dtype=np.float64)
            rd.read_dosages(int(i), d)
            d[d < 0] = np.nan
            cache[i] = d
        return cache[i]

    def r2(x, y, mask):
        x, y = x[mask], y[mask]
        ok = ~np.isnan(x) & ~np.isnan(y)
        if ok.sum() < 20 or x[ok].std() == 0 or y[ok].std() == 0:
            return np.nan
        return np.corrcoef(x[ok], y[ok])[0, 1] ** 2

    allm = np.ones(n, bool)
    rows = []
    for chrom, qc in q.assign(chrom=q["query"].str.split(":").str[0],
                              pos=q["query"].str.split(":").str[1].astype(int)).groupby("chrom"):
        pc = pv[pv.chrom == chrom]
        gc = gw[gw.chrom == chrom]
        pos2idx = pc.groupby("pos").idx.apply(list).to_dict()
        for x in qc.itertuples():
            if x.query not in id2idx:
                continue
            qd = dos(id2idx[x.query])
            near = gc[(gc.pos - x.pos).abs() <= WINDOW]
            for site, assoc in near.groupby("pos"):
                v = np.array([[r2(qd, dos(i), allm), r2(qd, dos(i), eur)] for i in pos2idx.get(site, [])])
                if not len(v) or np.isnan(v).all():
                    continue  # GWAS position not in the cohort panel
                best = np.nanmax(np.where(np.isnan(v), -1, v), axis=0)  # best allele at that position
                best[best < 0] = np.nan
                if np.nanmax(best) < 0.5:
                    continue
                for s in assoc.itertuples():
                    rows.append({"gene_ncbi": x.gene_ncbi, "role": x.role, "query": x.query, "gwas_pos": site,
                                 "r2_all": best[0], "r2_eur": best[1], "snp": s.snp, "risk_allele": s.risk_allele,
                                 "trait": s.trait, "mapped_trait": s.mapped_trait, "mlog10p": s.mlog10p,
                                 "or_beta": s.or_beta, "study": s.study, "pubmed": s.pubmed, "author": s.author})
        cache.clear()
    hits = pd.DataFrame(rows)
    hits = hits.merge(g[["gene_ncbi", "symbol", "tier", "region_class"]], on="gene_ncbi")
    hits["r2_max"] = hits[["r2_all", "r2_eur"]].max(axis=1)
    hits["proxy"] = hits.r2_max >= 0.8
    hits["genome_wide"] = hits.mlog10p >= -np.log10(5e-8)
    hits["immune_or_blood"] = hits.trait.fillna("").str.contains(IMMUNE) | hits.mapped_trait.fillna("").str.contains(IMMUNE)
    hits.sort_values(["tier", "symbol", "r2_max"], ascending=[True, True, False]).to_csv(
        out / f"lookup-{args.cohort}-{a}-vs-{b}{tag}.tsv", sep="\t", index=False, float_format="%.4g")

    sig = hits[hits.proxy & hits.genome_wide]
    per = []
    for r in g.itertuples():
        s = sig[sig.gene_ncbi == r.gene_ncbi]
        roles = lambda k: s[s.role == k]
        per.append({"gene_ncbi": r.gene_ncbi, "symbol": r.symbol, "tier": r.tier, "region_class": r.region_class,
                    "opposite_sign": np.sign(r.slope_a) != np.sign(r.slope_b),
                    "gwas_hits": s.study.nunique(), "traits": "; ".join(sorted(set(s.trait.dropna()))[:12]),
                    "immune_or_blood": bool(s.immune_or_blood.any()),
                    "via_pair": bool(len(s[s.role.str.startswith("pair")])),
                    f"via_lead_{a}": bool(len(roles(f"lead_{a}"))), f"via_lead_{b}": bool(len(roles(f"lead_{b}"))),
                    "best_r2": s.r2_max.max() if len(s) else np.nan,
                    "best_mlog10p": s.mlog10p.max() if len(s) else np.nan})
    per = pd.DataFrame(per)
    per.to_csv(out / f"genes-{args.cohort}-{a}-vs-{b}{tag}.tsv", sep="\t", index=False, float_format="%.4g")

    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 90)
    for lab, tiers in (("core", ["core"]), ("main", ["core", "main"]), ("extent", ["core", "main", "extent"]),
                       ("background", ["background"])):
        if not per.tier.isin(tiers).any():
            continue
        x = per[per.tier.isin(tiers)]
        print(f"{lab}: {len(x)} genes; with a genome-wide GWAS hit at r2 >= 0.8: {int((x.gwas_hits > 0).sum())}, "
              f"immune or blood: {int(x.immune_or_blood.sum())}; only via one reference's own lead: "
              f"{int(((x[f'via_lead_{a}'] != x[f'via_lead_{b}']) & (x.gwas_hits > 0)).sum())}")
    x = per[per.tier.isin(["core", "main"]) & (per.gwas_hits > 0)].sort_values("best_mlog10p", ascending=False)
    print(x[["symbol", "tier", "region_class", "opposite_sign", "gwas_hits", "immune_or_blood", "via_pair",
             f"via_lead_{a}", f"via_lead_{b}", "best_r2", "best_mlog10p", "traits"]].to_string(index=False))


if __name__ == "__main__":
    main()
