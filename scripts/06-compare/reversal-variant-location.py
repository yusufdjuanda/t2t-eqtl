"""Step 6: where do the lead variants of reversed eQTLs lie? (robust genes with opposite-sign effects at the same variant)

    pixi run python scripts/06-compare/reversal-variant-location.py --cohort mage [--a chm13v2.0] [--b grch38-ffx]

Why: a variant inside the gene's transcribed sequence can steer reads between near-identical copies during alignment
(a "mapping eQTL"), which would explain a sign flip between references without any regulatory effect. For every
robust gene (no artefact flag) whose slopes have opposite signs, the compared variant is located on each reference:
in an exon of the gene, in its body (intron), in an exon of its read partner (read-destination table), in exons of
other genes, the distance to the gene, and on T2T the segmental duplication (SD) overlap (>= 98% identity = UCSC
colours orange/yellow, as in make-gene-classes.py). Reversal class from the nominal p of each slope: "both strong"
(p < 1e-5 on both), "both nominal" (p < 0.01 on both), "one side ~zero".
Writes results/06-compare/eqtl-comparison/<cohort>/robust-matched/reversal-variant-location-<A>-vs-<B>.tsv.
"""
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SD98 = {"204,204,0", "255,103,0"}


def exons(ref, chroms):
    """Exon records of one reference as chrom/start/end/gene_ncbi (1-based closed), chromosomes in chrN style."""
    alias = {}
    if ref == "chm13v2.0":
        m = pd.read_csv(ROOT / "config/chrom-names-chm13v2.0.tsv", sep="\t")
        alias = dict(zip(m.refseq, m.ucsc))
    keep = {c for c, u in alias.items() if u in chroms} if alias else set(chroms)
    rows = []
    gid = re.compile(r'GeneID:(\d+)')
    with open(ROOT / f"data/star-references/{ref}/annotation.gtf") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split("\t", 8)
            if len(p) < 9 or p[2] != "exon" or p[0] not in keep:
                continue
            g = gid.search(p[8])
            if g:
                rows.append((alias.get(p[0], p[0]), int(p[3]), int(p[4]), g.group(1)))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "gene_ncbi"]).drop_duplicates()


def parse_variant(v):
    c, pos, ref, _ = v.split(":")
    return c, int(pos), int(pos) + len(ref) - 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--a", default="chm13v2.0")
    ap.add_argument("--b", default="grch38-ffx")
    args = ap.parse_args()
    out = ROOT / f"results/06-compare/eqtl-comparison/{args.cohort}/robust-matched"
    r = pd.read_csv(out / f"changed-{args.a}-vs-{args.b}.tsv", sep="\t", dtype={"gene_ncbi": str})
    r = r[~r.artefact_flag & (np.sign(r.slope_a) != np.sign(r.slope_b))].copy()
    r["reversal"] = np.select([(r.log10p_a >= 5) & (r.log10p_b >= 5), (r.log10p_a >= 2) & (r.log10p_b >= 2)],
                              ["both strong", "both nominal"], "one side ~zero")
    rd = pd.read_csv(out / f"read-destination-{args.a}-vs-{args.b}.tsv", sep="\t", dtype={"gene_ncbi": str})
    r = r.merge(rd[["gene_ncbi", "partner", "partner_r"]], on="gene_ncbi", how="left")

    res = r[["gene_ncbi", "symbol", "biotype", "region_class", "reversal", "slope_a", "slope_b", "log10p_a", "log10p_b",
             "z_diff", "z_diff_salmon", "partner", "partner_r"]].copy()
    for side, ref, vcol in [("a", args.a, "var_a"), ("b", args.b, "var_b")]:
        genes = pd.read_csv(ROOT / f"data/star-references/{ref}/genes.tsv", sep="\t", dtype={"gene_ncbi": str})
        if ref == "chm13v2.0":
            m = pd.read_csv(ROOT / "config/chrom-names-chm13v2.0.tsv", sep="\t")
            genes["chrom"] = genes.chrom.map(dict(zip(m.refseq, m.ucsc))).fillna(genes.chrom)
        sym2id = dict(zip(genes.symbol, genes.gene_ncbi))
        id2sym = dict(zip(genes.gene_ncbi, genes.symbol))
        var = r[vcol].map(parse_variant)
        ex = exons(ref, set(v[0] for v in var))
        loc = []
        for (c, s, e), g, partner in zip(var, r.gene_ncbi, r.partner):
            hit = ex[(ex.chrom == c) & (ex.start <= e) & (ex.end >= s)]
            body = genes[(genes.gene_ncbi == g) & (genes.chrom == c)]
            dist = np.nan if body.empty else int(max(0, body.start.min() - e, s - body.end.max()))
            pid = sym2id.get(partner)
            others = sorted(set(id2sym.get(x, x) for x in hit.gene_ncbi) - {id2sym.get(g)})
            loc.append({f"variant_{side}": f"{c}:{s}", f"in_gene_exon_{side}": g in set(hit.gene_ncbi),
                        f"in_gene_body_{side}": dist == 0, f"dist_to_gene_kb_{side}": round(dist / 1000, 1),
                        f"in_partner_exon_{side}": pid is not None and pid in set(hit.gene_ncbi),
                        f"other_exons_{side}": ",".join(others)})
        res = pd.concat([res, pd.DataFrame(loc, index=res.index)], axis=1)

    sd = pd.read_csv(ROOT / "data/references/t2t-annotations/chm13v2.0_SD.bed", sep="\t", header=None,
                     usecols=[0, 1, 2, 8], names=["chrom", "start", "end", "rgb"])
    sd_any, sd98 = [], []
    for c, s, e in r.var_a.map(parse_variant):
        h = sd[(sd.chrom == c) & (sd.start < e) & (sd.end >= s)]
        sd_any.append(not h.empty)
        sd98.append(bool(h.rgb.isin(SD98).any()))
    res["in_sd_a"], res["in_sd98_a"] = sd_any, sd98
    res = res.sort_values(["reversal", "z_diff"], key=lambda x: x if x.name == "reversal" else -x.abs())
    res.to_csv(out / f"reversal-variant-location-{args.a}-vs-{args.b}.tsv", sep="\t", index=False, float_format="%.4g")
    pd.set_option("display.width", 250)
    cols = ["symbol", "reversal", "partner", "in_gene_exon_a", "in_gene_exon_b", "in_gene_body_a", "dist_to_gene_kb_a",
            "in_partner_exon_a", "in_partner_exon_b", "in_sd_a", "in_sd98_a", "other_exons_a"]
    print(res[cols].to_string(index=False))
    print(res.groupby("reversal")[["in_gene_exon_a", "in_gene_exon_b", "in_gene_body_a", "in_partner_exon_a",
                                   "in_sd_a", "in_sd98_a"]].sum())


if __name__ == "__main__":
    main()
