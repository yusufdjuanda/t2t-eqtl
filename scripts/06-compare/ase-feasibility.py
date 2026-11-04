"""Step 6: can allele-specific expression (ASE) decide the direction of the changed eQTLs? Feasibility, no new data.

    pixi run python scripts/06-compare/ase-feasibility.py [--cohort mage] [--threads 8]   (needs minimap2 on PATH)

ASE test of an eQTL: in a person heterozygous at the lead variant, the haplotype carrying the alt allele should give
more (slope > 0) or fewer (slope < 0) RNA reads than the other haplotype. The reads are told apart at a heterozygous
SNP in the gene's exons, phased with the lead in the panel. This only works if a read covering that exon SNP can be
placed on the gene and not on a paralog, which is the problem for most of the changed genes.
For the main gene set (T2T vs patched GRCh38, tiers core and main of gene-sets/changed-*.tsv, define-gene-sets.py;
117 genes) this counts, per gene:
  lead_het                  individuals heterozygous at the T2T lead variant
  exon_snps                 biallelic SNPs of the T2T panel in the gene's T2T exons
  distinct_t2t / _ffx       of these, SNPs whose 150 bp window (both alleles) maps to its own position, with the best
                            other hit at least 2 mismatches worse (minimap2 -x sr, alignment score gap >= 20;
                            "_1mm": gap >= 10), on T2T and on patched GRCh38 (position via the variant matching table)
  informative_t2t / _ffx / _both  individuals heterozygous at the lead and at >= 1 distinct exon SNP
  mean_count_t2t / _ffx, depth_t2t / _ffx  mean STAR count and approximate reads per exonic base per individual
Writes results/06-compare/ase/feasibility-<cohort>.tsv and the per-SNP table snps-<cohort>.tsv.gz.
"""
import argparse
import os
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pgenlib

ROOT = Path(__file__).resolve().parents[2]
HALF = 75  # window = 150 bp, the MAGE read length


def exons(gtf, genes, rename=None):
    """Exon intervals (chrom, start, end; 1-based inclusive) of the given GeneIDs."""
    rows = []
    with open(gtf) as f:
        for line in f:
            if "\texon\t" not in line:
                continue
            m = re.search(r'db_xref "GeneID:(\d+)"', line)
            if m and m.group(1) in genes:
                c = line.split("\t")
                rows.append((m.group(1), rename.get(c[0], c[0]) if rename else c[0], int(c[3]), int(c[4])))
    return pd.DataFrame(rows, columns=["gene_ncbi", "chrom", "start", "end"]).drop_duplicates()


def reader(fa):
    """Function (chrom, start, end) -> sequence (0-based, half-open) from an indexed FASTA."""
    fai = {}
    for l in open(f"{fa}.fai"):
        c = l.split("\t")
        fai[c[0]] = (int(c[1]), int(c[2]), int(c[3]), int(c[4]))
    fh = open(fa, "rb")

    def get(chrom, start, end):
        length, off, lb, lw = fai[chrom]
        start, end = max(start, 0), min(end, length)
        a = off + start // lb * lw + start % lb
        b = off + (end - 1) // lb * lw + (end - 1) % lb + 1
        fh.seek(a)
        return fh.read(b - a).replace(b"\n", b"").replace(b"\r", b"").decode().upper()
    return get


def windows(snps, get, tag):
    """FASTA records of the 150 bp windows around each SNP, ref and alt allele."""
    out, bad = [], 0
    for r in snps.itertuples():
        s = get(r.chrom, r.pos - 1 - HALF, r.pos + HALF - 1)
        if len(s) != 2 * HALF or s[HALF] != r.ref:
            bad += 1
            continue
        for lab, al in (("ref", r.ref), ("alt", r.alt)):
            out.append(f">{tag}|{r.snp}|{lab}|{r.chrom}|{r.pos}\n{s[:HALF]}{al}{s[HALF + 1:]}\n")
    return "".join(out), bad


def distinct(fasta, genome, threads, tmp):
    """Per SNP: maps to its own position with both alleles, and the score gap to the best other hit."""
    q = tmp / "windows.fa"
    q.write_text(fasta)
    sam = subprocess.run(["minimap2", "-t", str(threads), "-x", "sr", "-a", "--secondary=yes", "-N", "50", "-p", "0.3",
                          "-I", "8G", genome, str(q)], check=True, capture_output=True, text=True).stdout
    best = {}
    for line in sam.splitlines():
        if line.startswith("@"):
            continue
        c = line.split("\t")
        flag = int(c[1])
        if flag & 4 or flag & 2048:
            continue
        name = c[0]
        tag, snp, lab, chrom, pos = name.split("|")
        a = int(re.search(r"AS:i:(-?\d+)", line).group(1))
        here = c[2] == chrom and abs(int(c[3]) - (int(pos) - HALF)) <= 10
        d = best.setdefault(name, {"own": None, "other": None})
        if here:
            d["own"] = a if d["own"] is None else max(d["own"], a)
        else:
            d["other"] = a if d["other"] is None else max(d["other"], a)
    rows = []
    for name, d in best.items():
        tag, snp, lab, chrom, pos = name.split("|")
        gap = np.nan if d["own"] is None else (d["own"] - d["other"] if d["other"] is not None else 300)
        rows.append((snp, lab, gap))
    g = pd.DataFrame(rows, columns=["snp", "allele", "gap"]).pivot(index="snp", columns="allele", values="gap")
    return g.min(axis=1)  # the worse allele decides


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", default="mage")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()
    out = ROOT / "results/06-compare/ase"
    out.mkdir(parents=True, exist_ok=True)
    tmp = Path(os.environ["TMPDIR"]) / "ase-feasibility"
    tmp.mkdir(parents=True, exist_ok=True)

    base = ROOT / f"results/06-compare/eqtl-comparison/{args.cohort}/robust-matched"
    g = pd.read_csv(ROOT / f"results/06-compare/eqtl-comparison/{args.cohort}/gene-sets/changed-chm13v2.0-vs-grch38-ffx.tsv",
                    sep="\t", dtype={"gene_ncbi": str})
    g = g[g.tier.isin(["core", "main"])]
    rev = pd.read_csv(base / "reversal-variant-location-chm13v2.0-vs-grch38-ffx.tsv", sep="\t", dtype={"gene_ncbi": str})
    g["true_reversal"] = g.gene_ncbi.isin(rev.loc[rev.reversal.str.contains("strong"), "gene_ncbi"])
    genes = set(g.gene_ncbi)
    print(f"{len(g)} genes, {int(g.true_reversal.sum())} true reversals", flush=True)

    names = pd.read_csv(ROOT / "config/chrom-names-chm13v2.0.tsv", sep="\t")
    nc2chr = dict(zip(names.refseq, names.ucsc)); chr2nc = dict(zip(names.ucsc, names.refseq))
    ex_t = exons(ROOT / "data/star-references/chm13v2.0/annotation.gtf", genes, nc2chr)
    ex_f = exons(ROOT / "data/star-references/grch38-ffx/annotation.gtf", genes)

    # T2T panel SNPs in the exons
    gdir = ROOT / f"results/05-eqtl/{args.cohort}/genotypes"
    chroms = set(ex_t.chrom)
    pv = pd.read_csv(gdir / "t2t.all.pvar", sep="\t", usecols=[0, 1, 2, 3, 4], names=["chrom", "pos", "id", "ref", "alt"],
                     comment="#", dtype={"chrom": str})
    pv["idx"] = np.arange(len(pv))
    lead_idx = dict(zip(pv.id, pv.idx))
    pv = pv[pv.chrom.isin(chroms) & (pv.ref.str.len() == 1) & (pv.alt.str.len() == 1)]
    hits = []
    for e in ex_t.itertuples():
        x = pv[(pv.chrom == e.chrom) & (pv.pos >= e.start) & (pv.pos <= e.end)]
        hits.append(x.assign(gene_ncbi=e.gene_ncbi))
    snps = pd.concat(hits).drop_duplicates(["gene_ncbi", "id"]).rename(columns={"id": "snp"})
    print(f"{len(snps)} exon SNPs in {snps.gene_ncbi.nunique()} genes", flush=True)

    # GRCh38 position of each SNP (matched variants only)
    m = pd.concat(pd.read_csv(f, sep="\t", usecols=["side", "category", "t2t_id", "grch38_id"], low_memory=False)
                  for f in sorted((ROOT / "results/06-compare/variant-matching").glob("chr*.tsv.gz")))
    m = m[(m.side == "t2t") & (m.category == "matched")].set_index("t2t_id").grch38_id
    snps["grch38_id"] = snps.snp.map(m)

    # distinctness on each reference
    get_t = reader(ROOT / "data/star-references/chm13v2.0/genome.fa")
    fa, bad = windows(snps.assign(chrom=snps.chrom.map(chr2nc)), get_t, "t2t")
    print(f"T2T windows: {len(snps) - bad} (ref base mismatch {bad})", flush=True)
    snps["gap_t2t"] = snps.snp.map(distinct(fa, str(ROOT / "data/star-references/chm13v2.0/genome.fa"), args.threads, tmp))
    f = snps.dropna(subset=["grch38_id"]).copy()
    p = f.grch38_id.str.split(":")
    f = f.assign(chrom=p.str[0], pos=p.str[1].astype(int), ref=p.str[2], alt=p.str[3], snp=f.snp)
    f = f[(f.ref.str.len() == 1) & (f.alt.str.len() == 1)]
    get_f = reader(ROOT / "data/star-references/grch38-ffx/genome.fa")
    fa, bad = windows(f, get_f, "ffx")
    print(f"GRCh38-ffx windows: {len(f) - bad} (ref base mismatch {bad})", flush=True)
    snps["gap_ffx"] = snps.snp.map(distinct(fa, str(ROOT / "data/star-references/grch38-ffx/genome.fa"), args.threads, tmp))
    in_ffx_exon = []
    for r in snps.itertuples():
        if not isinstance(r.grch38_id, str):
            in_ffx_exon.append(False)
            continue
        c, pos = r.grch38_id.split(":")[0], int(r.grch38_id.split(":")[1])
        e = ex_f[(ex_f.gene_ncbi == r.gene_ncbi) & (ex_f.chrom == c)]
        in_ffx_exon.append(bool(((e.start <= pos) & (e.end >= pos)).any()))
    snps["in_ffx_exon"] = in_ffx_exon
    for s in ("t2t", "ffx"):
        snps[f"distinct_{s}"] = snps[f"gap_{s}"] >= 20
        snps[f"distinct_{s}_1mm"] = snps[f"gap_{s}"] >= 10
    snps["distinct_ffx"] &= snps.in_ffx_exon
    snps["distinct_ffx_1mm"] &= snps.in_ffx_exon

    # phased genotypes
    rd = pgenlib.PgenReader(str(gdir / "t2t.all.pgen").encode())
    n = rd.get_raw_sample_ct()
    psam = pd.read_csv(gdir / "t2t.all.psam", sep="\t")
    def het(i):
        a = np.empty(2 * n, dtype=np.int32); ph = np.empty(n, dtype=np.uint8)
        rd.read_alleles_and_phasepresent(int(i), a, ph)
        return a[0::2] != a[1::2]

    # expression
    runs = pd.read_csv(ROOT / f"results/04-expression/{args.cohort}/matrices/runs.tsv", sep="\t")
    gt = pd.read_csv(ROOT / "data/star-references/chm13v2.0/genes.tsv", sep="\t", dtype={"gene_ncbi": str})
    gf = pd.read_csv(ROOT / "data/star-references/grch38-ffx/genes.tsv", sep="\t", dtype={"gene_ncbi": str})
    cnt = {}
    for ref, gtab in (("t2t", gt), ("ffx", gf)):
        r = {"t2t": "chm13v2.0", "ffx": "grch38-ffx"}[ref]
        c = pd.read_csv(ROOT / f"results/04-expression/{args.cohort}/matrices/{r}/star-counts.tsv.gz", sep="\t",
                        index_col=0)
        ids = gtab.set_index("gene_id").gene_ncbi
        c = c[c.index.isin(ids.index[ids.isin(genes)])]
        c.index = ids.reindex(c.index).values
        cnt[ref] = c.mean(axis=1)
    rl = runs.read_length.median()

    def exlen(ex, gene):
        e = ex[ex.gene_ncbi == gene].sort_values(["chrom", "start"])
        tot, cur = 0, None
        for r in e.itertuples():
            if cur is None or r.chrom != cur[0] or r.start > cur[2]:
                if cur: tot += cur[2] - cur[1] + 1
                cur = [r.chrom, r.start, r.end]
            else:
                cur[2] = max(cur[2], r.end)
        return tot + (cur[2] - cur[1] + 1 if cur else 0)

    rows = []
    for r in g.itertuples():
        rec = {"gene_ncbi": r.gene_ncbi, "symbol": r.symbol, "region_class": r.region_class,
               "tier": r.tier, "true_reversal": r.true_reversal,
               "opposite_sign": np.sign(r.slope_a) != np.sign(r.slope_b), "slope_t2t": r.slope_a, "slope_ffx": r.slope_b, "lead": r.var_a}
        if r.var_a not in lead_idx:
            rows.append(rec); continue
        lh = het(lead_idx[r.var_a])
        s = snps[snps.gene_ncbi == r.gene_ncbi]
        H = {x.snp: het(x.idx) for x in s.itertuples()}
        rec.update(lead_het=int(lh.sum()), exon_snps=len(s))
        for k in ("t2t", "ffx"):
            for suf in ("", "_1mm"):
                d = s[s[f"distinct_{k}{suf}"]]
                rec[f"distinct_{k}{suf}"] = len(d)
                any_het = np.zeros(n, bool)
                for x in d.snp:
                    any_het |= H[x]
                rec[f"informative_{k}{suf}"] = int((lh & any_het).sum())
        both = s[s.distinct_t2t & s.distinct_ffx]
        any_het = np.zeros(n, bool)
        for x in both.snp:
            any_het |= H[x]
        rec["informative_both"] = int((lh & any_het).sum())
        for k, ex in (("t2t", ex_t), ("ffx", ex_f)):
            L = exlen(ex, r.gene_ncbi)
            rec[f"mean_count_{k}"] = cnt[k].get(r.gene_ncbi, np.nan)
            rec[f"depth_{k}"] = rec[f"mean_count_{k}"] * 2 * rl / L if L else np.nan
        rows.append(rec)
    t = pd.DataFrame(rows)
    t.to_csv(out / f"feasibility-{args.cohort}.tsv", sep="\t", index=False, float_format="%.4g")
    snps.drop(columns=["idx"]).to_csv(out / f"snps-{args.cohort}.tsv.gz", sep="\t", index=False, float_format="%.4g")

    pd.set_option("display.width", 250)
    cols = ["symbol", "tier", "true_reversal", "opposite_sign", "region_class", "lead_het", "exon_snps", "distinct_t2t", "distinct_ffx",
            "informative_t2t", "informative_ffx", "informative_both", "informative_t2t_1mm", "depth_t2t", "depth_ffx"]
    print(t.sort_values(["true_reversal", "informative_both"], ascending=False)[cols].to_string(index=False))
    for lab, x in (("true reversals", t[t.true_reversal]), ("core", t[t.tier == "core"]), ("main (incl. core)", t)):
        print(f"\n{lab}: {len(x)} genes; >= 20 informative individuals on both references: "
              f"{int((x.informative_both >= 20).sum())}; >= 10: {int((x.informative_both >= 10).sum())}; "
              f"on T2T only (>= 20): {int((x.informative_t2t >= 20).sum())}")


if __name__ == "__main__":
    main()
