"""Step 6: where do the reads of a changed gene go on the other reference? (run after combine-methods.py)

    pixi run python scripts/06-compare/read-destination.py --cohort mage [--a chm13v2.0] [--b grch38-ffx] [--method star]

For each robustly changed gene (results/06-compare/eqtl-comparison/<cohort>/robust-matched/changed-<A>-vs-<B>.tsv,
no artefact flag), per individual: delta = CPM on A - CPM on B (genes matched by NCBI GeneID; a gene missing from
one reference's matrix counts as 0 there, so reads gained by an annotation-only gene are seen too). Reads that one
reference gives gene X and the other gives gene Y make delta_X and delta_Y move in opposite directions across
individuals: the partner is the gene whose delta correlates most negatively with delta_X (|mean delta_Y| > 0.1
CPM). Reported per gene: mean CPM on A and B, mean delta, the best partner (symbol, class, correlation, and the
fraction of X's per-individual delta it accounts for = -slope of delta_Y on delta_X), the next two partners.
Writes results/06-compare/eqtl-comparison/<cohort>/robust-matched/read-destination-<A>-vs-<B>.tsv.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def cpm_by_geneid(cohort, ref, method):
    c = pd.read_csv(ROOT / f"results/03-expression/{cohort}/matrices/{ref}/{method}-counts.tsv.gz", sep="\t", index_col=0)
    g = pd.read_csv(ROOT / f"data/star-references/{ref}/genes.tsv", sep="\t", dtype={"gene_ncbi": str})
    c = c.groupby(c.index.map(dict(zip(g.gene_id, g.gene_ncbi)))).sum()
    return c / c.sum() * 1e6


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--a", default="chm13v2.0")
    ap.add_argument("--b", default="grch38-ffx")
    ap.add_argument("--method", default="star")
    args = ap.parse_args()
    out = ROOT / f"results/06-compare/eqtl-comparison/{args.cohort}/robust-matched"
    r = pd.read_csv(out / f"changed-{args.a}-vs-{args.b}.tsv", sep="\t", dtype={"gene_ncbi": str})
    r = r[~r.artefact_flag]
    a, b = cpm_by_geneid(args.cohort, args.a, args.method), cpm_by_geneid(args.cohort, args.b, args.method)
    ids = a.columns.intersection(b.columns)
    genes = a.index.union(b.index)
    a, b = a.reindex(index=genes, columns=ids, fill_value=0), b.reindex(index=genes, columns=ids, fill_value=0)
    d = a - b
    d = d[(d.mean(axis=1).abs() > 0.1) | d.index.isin(r.gene_ncbi)]
    cls = pd.read_csv(ROOT / "results/06-compare/gene-classes.tsv", sep="\t", dtype={"gene_ncbi": str}).set_index("gene_ncbi")
    dc = d.sub(d.mean(axis=1), axis=0)
    dn = dc.div(np.sqrt((dc ** 2).sum(axis=1)), axis=0).fillna(0)
    rows = []
    for _, x in r.iterrows():
        g = x.gene_ncbi
        if g not in d.index:
            continue
        cor = dn.drop(index=g) @ dn.loc[g]
        best = cor.nsmallest(3)
        v = dc.loc[g]; vv = (v ** 2).sum()
        frac = {p: float(-(dc.loc[p] @ v) / vv) if vv > 0 else np.nan for p in best.index}
        sym = lambda p: cls.symbol.get(p, p)
        rows.append({"gene_ncbi": g, "symbol": x.symbol, "region_class": x.region_class, "stronger_on": x.stronger_on,
                     f"cpm_{args.a}": a.loc[g].mean(), f"cpm_{args.b}": b.loc[g].mean(), "mean_delta": d.loc[g].mean(),
                     "partner": sym(best.index[0]), "partner_class": cls.region_class.get(best.index[0], "NA"),
                     "partner_annotated_on": cls.annotated_on.get(best.index[0], "NA"),
                     "partner_r": best.iloc[0], "partner_frac": frac[best.index[0]],
                     "next_partners": "; ".join(f"{sym(p)} (r={c:.2f})" for p, c in best.iloc[1:].items())})
    t = pd.DataFrame(rows).sort_values("partner_r")
    t.to_csv(out / f"read-destination-{args.a}-vs-{args.b}.tsv", sep="\t", index=False, float_format="%.4g")
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
    print(t.drop(columns=["gene_ncbi"]).round(2).to_string(index=False))
    print(f"\n{len(t)} genes; partner with r < -0.5: {(t.partner_r < -0.5).sum()}; r < -0.8: {(t.partner_r < -0.8).sum()}")


if __name__ == "__main__":
    main()
