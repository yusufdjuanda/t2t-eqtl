"""Step 6: label every gene (NCBI GeneID) with the region class used to test H0, plus artefact flags.

    pixi run python scripts/06-compare/make-gene-classes.py

Inputs: data/star-references/{chm13v2.0,grch38}/genes.tsv (grch38-ffx uses the grch38 genes); FixItFelix Additional
file 2 Tables S1 (falsely duplicated: correct and false copies) and S2 (collapsed), GRCh38 coordinates;
data/references/t2t-annotations/chm13v2.0_SD.bed; data/external/t2t-consortium/chm13v2-unique_to_hg38.bed
(T2T sequence absent from GRCh38; UCSC names, mapped with config/chrom-names-chm13v2.0.tsv).
Writes results/06-compare/gene-classes.tsv, one row per GeneID:
  annotated_on      both | t2t_only | grch38_only
  ffx_patched       a GRCh38 locus overlaps a FixItFelix patched region (either copy of a false duplication)
  t2t_sd_frac       largest fraction of a T2T locus inside segmental duplications of any identity (>= 90%)
  t2t_sd98_frac     the same for high-identity SDs only (>= 98%: UCSC colours yellow 204,204,0 = 98-99% and orange
                    255,103,0 = > 99%; grey shades = 90-98%). Grey SDs include old copies such as processed-pseudogene
                    parents (ACTB is 73% covered by grey SDs only), which short reads separate easily.
  t2t_new_frac      largest fraction of a T2T locus inside T2T sequence absent from GRCh38
  region_class      first that applies: ffx_patched > t2t_new_sequence (t2t_new_frac >= 0.5) > sd (t2t_sd98_frac >= 0.5) > other
  flag_mhc, flag_par, flag_chrM, flag_rrna    genes to report separately (MHC haplotype differs between references;
                    PAR copies; mitochondrial; rRNA). Annotation-difference artefacts (e.g. HCP5, NEAT1, P2RY8) are
                    flagged after the expression data are in (step 6), not here.
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def merged(iv):
    """dict chrom -> (starts, ends, cumulative lengths) of merged half-open intervals."""
    out = {}
    for c, d in iv.groupby("chrom"):
        m = []
        for s, e in sorted(zip(d.start, d.end)):
            if m and s <= m[-1][1]:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        a = np.array(m)
        out[c] = (a[:, 0], a[:, 1], np.concatenate([[0], np.cumsum(a[:, 1] - a[:, 0])]))
    return out


def covered_frac(idx, chrom, s, e):
    """fraction of [s, e) covered by the merged intervals on chrom."""
    if chrom not in idx or e <= s:
        return 0.0
    st, en, cum = idx[chrom]

    def upto(x):
        i = np.searchsorted(st, x, side="right")
        return 0 if i == 0 else cum[i - 1] + min(max(x - st[i - 1], 0), en[i - 1] - st[i - 1])
    return (upto(e) - upto(s)) / (e - s)


def main():
    names = pd.read_csv(ROOT / "config/chrom-names-chm13v2.0.tsv", sep="\t")
    to_refseq = dict(zip(names.ucsc, names.refseq))

    def ucsc_bed(path, colours=None):
        b = pd.read_csv(path, sep="\t", header=None, comment="#")
        if colours is not None:
            b = b[b[8].isin(colours)]
        b = b.iloc[:, :3]
        b.columns = ["chrom", "start", "end"]
        b["chrom"] = b.chrom.map(to_refseq)
        return merged(b.dropna())
    sd = ucsc_bed(ROOT / "data/references/t2t-annotations/chm13v2.0_SD.bed")
    sd98 = ucsc_bed(ROOT / "data/references/t2t-annotations/chm13v2.0_SD.bed", colours={"204,204,0", "255,103,0"})
    new = ucsc_bed(ROOT / "data/external/t2t-consortium/chm13v2-unique_to_hg38.bed")

    x = pd.ExcelFile(ROOT / "data/external/fixitfelix-behera2023/13059_2023_2863_MOESM2_ESM.xlsx")
    iv = []
    for v in x.parse("Table S1 Duplicated regions", header=None).iloc[3:, :2].values.ravel():
        if isinstance(v, str) and re.match(r"chr\S+:\s*\d", v):
            c, r = v.split(":")
            a, b = [int(t.strip().replace(",", "")) for t in r.split("-")]
            iv.append((c.strip(), a - 1, b))
    for _, r in x.parse("Table S2 Collapsed regions", header=None).iloc[3:12].iterrows():
        if isinstance(r[0], str) and r[0].startswith("chr"):
            iv.append((r[0], int(r[1]) - 1, int(r[2])))
    ffx = merged(pd.DataFrame(iv, columns=["chrom", "start", "end"]))

    t2t = pd.read_csv(ROOT / "data/star-references/chm13v2.0/genes.tsv", sep="\t", dtype={"gene_ncbi": str})
    g38 = pd.read_csv(ROOT / "data/star-references/grch38/genes.tsv", sep="\t", dtype={"gene_ncbi": str})
    t2t["sd"] = [covered_frac(sd, c, s - 1, e) for c, s, e in zip(t2t.chrom, t2t.start, t2t.end)]
    t2t["sd98"] = [covered_frac(sd98, c, s - 1, e) for c, s, e in zip(t2t.chrom, t2t.start, t2t.end)]
    t2t["new"] = [covered_frac(new, c, s - 1, e) for c, s, e in zip(t2t.chrom, t2t.start, t2t.end)]
    g38["ffx"] = [covered_frac(ffx, c, s - 1, e) > 0 for c, s, e in zip(g38.chrom, g38.start, g38.end)]

    def locus(d):
        return d.chrom + ":" + d.start.astype(str) + "-" + d.end.astype(str)
    a = t2t.assign(locus=locus(t2t)).groupby("gene_ncbi").agg(
        symbol_t2t=("symbol", "first"), biotype_t2t=("biotype", "first"), loci_t2t=("locus", ";".join),
        n_loci_t2t=("locus", "size"), t2t_sd_frac=("sd", "max"), t2t_sd98_frac=("sd98", "max"), t2t_new_frac=("new", "max"))
    b = g38.assign(locus=locus(g38)).groupby("gene_ncbi").agg(
        symbol_grch38=("symbol", "first"), biotype_grch38=("biotype", "first"), loci_grch38=("locus", ";".join),
        n_loci_grch38=("locus", "size"), ffx_patched=("ffx", "any"))
    t = a.join(b, how="outer")
    t["symbol"] = t.symbol_t2t.fillna(t.symbol_grch38)
    t["biotype"] = t.biotype_t2t.fillna(t.biotype_grch38)
    t["annotated_on"] = np.select([t.n_loci_t2t.notna() & t.n_loci_grch38.notna(), t.n_loci_t2t.notna()],
                                  ["both", "t2t_only"], "grch38_only")
    t["ffx_patched"] = t.ffx_patched.fillna(False).astype(bool)
    t[["t2t_sd_frac", "t2t_sd98_frac", "t2t_new_frac"]] = t[["t2t_sd_frac", "t2t_sd98_frac", "t2t_new_frac"]].fillna(0)
    t["region_class"] = np.select([t.ffx_patched, t.t2t_new_frac >= 0.5, t.t2t_sd98_frac >= 0.5],
                                  ["ffx_patched", "t2t_new_sequence", "sd"], "other")

    # flags: MHC span per reference (first to last HLA gene on the HLA chromosome), PARs, chrM, rRNA
    def mhc_span(d):
        h = d[d.symbol.str.match(r"^HLA-[A-Z]")]
        h = h[h.chrom == h.chrom.mode()[0]]
        return h.chrom.iloc[0], h.start.min(), h.end.max()
    mhc = set()
    for d in (t2t, g38):
        c, s, e = mhc_span(d)
        mhc |= set(d[(d.chrom == c) & (d.end >= s) & (d.start <= e)].gene_ncbi)
    par = set(g38[(g38.chrom.isin(["chrX", "chrY"])) & ((g38.end <= 2781479) | (g38.start >= 155701383))].gene_ncbi) | \
        set(g38[(g38.chrom == "chrY") & (g38.start >= 56887903)].gene_ncbi)
    t["flag_mhc"] = t.index.isin(mhc)
    t["flag_par"] = t.index.isin(par)
    t["flag_chrM"] = t.loci_grch38.fillna("").str.startswith("chrM:")
    t["flag_rrna"] = t.biotype == "rRNA"

    out = ROOT / "results/06-compare"
    out.mkdir(parents=True, exist_ok=True)
    cols = ["symbol", "biotype", "annotated_on", "region_class", "ffx_patched", "t2t_sd_frac", "t2t_sd98_frac", "t2t_new_frac",
            "flag_mhc", "flag_par", "flag_chrM", "flag_rrna", "n_loci_t2t", "n_loci_grch38", "loci_t2t", "loci_grch38"]
    t.index.name = "gene_ncbi"
    t[cols].to_csv(out / "gene-classes.tsv", sep="\t", float_format="%.3f")
    print(pd.crosstab(t.region_class, t.annotated_on, margins=True).to_string())
    pc = t.biotype == "protein_coding"
    print("\nprotein-coding:"); print(pd.crosstab(t[pc].region_class, t[pc].annotated_on, margins=True).to_string())
    print("\nflags:", {f: int(t[f].sum()) for f in ["flag_mhc", "flag_par", "flag_chrM", "flag_rrna"]})


if __name__ == "__main__":
    main()
