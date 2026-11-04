"""Step 6: targets for allele-specific expression (ASE): genes, regions and individuals.

    pixi run python scripts/06-compare/ase-targets.py [--cohort mage] [--min-informative 30] [--max-extra 120]

Genes: the main set (tiers core and main of gene-sets/changed-chm13v2.0-vs-grch38-ffx.tsv) and their read partners
(read-destination and reversal tables), so that reads moving between paralogs can be followed.
Regions:
  rna-<ref>.bed   gene bodies +- 5 kb on chm13v2.0 (RefSeq names, as the STAR index) and grch38-ffx; the RNA-seq
                  realignment keeps only reads here
  ont-t2t.bed     main genes +- 20 kb, and their T2T lead variant +- 20 kb, in the 1KGP ONT reference's names
                  (chr1...); long-read SNP calling and read-backed phasing there
Individuals: every cohort individual with a 1KGP ONT T2T alignment (long-read SNPs for the SD genes need their RNA as
well), plus up to --max-extra more, added greedily: each step takes the individual informative for the most
testable genes (ase-feasibility.py: heterozygous at the lead and at a distinguishable exon SNP) that still have
fewer than --min-informative people.
Writes results/06-compare/ase/targets/{genes.tsv, rna-chm13v2.0.bed, rna-grch38-ffx.bed, ont-t2t.bed,
individuals.tsv, runs.txt} and prints the expected informative counts.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pgenlib

ROOT = Path(__file__).resolve().parents[2]
ONT_LIST = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/t2t/"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", default="mage")
    ap.add_argument("--min-informative", type=int, default=30)
    ap.add_argument("--max-extra", type=int, default=120)
    ap.add_argument("--ont-samples", required=True, help="file with the IDs that have <id>.t2t.cram in " + ONT_LIST)
    args = ap.parse_args()
    out = ROOT / "results/06-compare/ase/targets"
    out.mkdir(parents=True, exist_ok=True)
    cmp_ = ROOT / f"results/06-compare/eqtl-comparison/{args.cohort}"

    g = pd.read_csv(cmp_ / "gene-sets/changed-chm13v2.0-vs-grch38-ffx.tsv", sep="\t", dtype={"gene_ncbi": str})
    g = g[g.tier.isin(["core", "main"])]
    rd = pd.read_csv(cmp_ / "robust-matched/read-destination-chm13v2.0-vs-grch38-ffx.tsv", sep="\t",
                     dtype={"gene_ncbi": str})
    rv = pd.read_csv(cmp_ / "robust-matched/reversal-variant-location-chm13v2.0-vs-grch38-ffx.tsv", sep="\t",
                     dtype={"gene_ncbi": str})
    partners = set(rd.partner.dropna()) | set(rv.partner.dropna())
    gt = pd.read_csv(ROOT / "data/star-references/chm13v2.0/genes.tsv", sep="\t", dtype={"gene_ncbi": str})
    gf = pd.read_csv(ROOT / "data/star-references/grch38-ffx/genes.tsv", sep="\t", dtype={"gene_ncbi": str})
    pid = set(gt.loc[gt.symbol.isin(partners), "gene_ncbi"]) | set(gf.loc[gf.symbol.isin(partners), "gene_ncbi"])
    genes = pd.DataFrame({"gene_ncbi": sorted(set(g.gene_ncbi) | pid)})
    genes["role"] = np.where(genes.gene_ncbi.isin(g.gene_ncbi), "target", "partner")
    sym = pd.concat([gt, gf]).drop_duplicates("gene_ncbi").set_index("gene_ncbi").symbol
    genes["symbol"] = genes.gene_ncbi.map(sym)
    genes.to_csv(out / "genes.tsv", sep="\t", index=False)
    print(f"{(genes.role == 'target').sum()} target genes, {(genes.role == 'partner').sum()} partners")

    for ref, tab in (("chm13v2.0", gt), ("grch38-ffx", gf)):
        b = tab[tab.gene_ncbi.isin(genes.gene_ncbi)][["chrom", "start", "end"]].copy()
        b["start"] = (b.start - 1 - 5000).clip(lower=0)
        b["end"] = b.end + 5000
        b.sort_values(["chrom", "start"]).to_csv(out / f"rna-{ref}.bed", sep="\t", header=False, index=False)

    names = pd.read_csv(ROOT / "config/chrom-names-chm13v2.0.tsv", sep="\t")
    nc2chr = dict(zip(names.refseq, names.ucsc))
    t = gt[gt.gene_ncbi.isin(g.gene_ncbi)][["chrom", "start", "end"]].copy()
    t["chrom"] = t.chrom.map(nc2chr)
    t["start"] = (t.start - 1 - 20000).clip(lower=0)
    t["end"] = t.end + 20000
    lead = g.var_a.str.split(":", expand=True)
    lv = pd.DataFrame({"chrom": lead[0], "start": (lead[1].astype(int) - 20000).clip(lower=0),
                       "end": lead[1].astype(int) + 20000})
    o = pd.concat([t, lv]).sort_values(["chrom", "start"])
    merged = []
    for c, x in o.groupby("chrom", sort=True):
        cur = None
        for r in x.itertuples():
            if cur and r.start <= cur[2]:
                cur[2] = max(cur[2], r.end)
            else:
                if cur: merged.append(cur)
                cur = [c, r.start, r.end]
        merged.append(cur)
    ont = pd.DataFrame(merged, columns=["chrom", "start", "end"])
    ont.to_csv(out / "ont-t2t.bed", sep="\t", header=False, index=False)
    print(f"ONT regions: {len(ont)}, {(ont.end - ont.start).sum() / 1e6:.1f} Mb")

    # individuals
    runs = pd.read_csv(ROOT / f"config/samplesheets/{args.cohort}.csv")
    cohort = set(runs["sample"])
    ont_ids = set(open(args.ont_samples).read().split()) & cohort
    snps = pd.read_csv(ROOT / f"results/06-compare/ase/snps-{args.cohort}.tsv.gz", sep="\t", dtype={"gene_ncbi": str})
    feas = pd.read_csv(ROOT / f"results/06-compare/ase/feasibility-{args.cohort}.tsv", sep="\t",
                       dtype={"gene_ncbi": str})
    test = feas[feas.informative_both >= 10]
    gdir = ROOT / f"results/05-eqtl/{args.cohort}/genotypes"
    pv = pd.read_csv(gdir / "t2t.all.pvar", sep="\t", comment="#", header=None, usecols=[2], names=["id"])
    idx = dict(zip(pv.id, range(len(pv))))
    psam = pd.read_csv(gdir / "t2t.all.psam", sep="\t").iloc[:, 0].values
    rdr = pgenlib.PgenReader(str(gdir / "t2t.all.pgen").encode())
    n = rdr.get_raw_sample_ct()

    def het(v):
        a = np.empty(2 * n, dtype=np.int32); p = np.empty(n, dtype=np.uint8)
        rdr.read_alleles_and_phasepresent(idx[v], a, p)
        return a[0::2] != a[1::2]

    inf = {}
    for r in test.itertuples():
        s = snps[(snps.gene_ncbi == r.gene_ncbi) & snps.distinct_t2t & snps.distinct_ffx]
        anyh = np.zeros(n, bool)
        for v in s.snp:
            anyh |= het(v)
        inf[r.symbol] = het(r.lead) & anyh
    M = pd.DataFrame(inf, index=psam)
    chosen = set(ont_ids)
    extra = []
    for _ in range(args.max_extra):
        have = M.loc[sorted(chosen)].sum()
        need = have[have < args.min_informative].index
        if not len(need):
            break
        rest = M.loc[~M.index.isin(chosen), need]
        gain = rest.sum(axis=1)
        if gain.max() == 0:
            break
        best = gain.idxmax()
        chosen.add(best); extra.append(best)
    ind = pd.DataFrame({"sample": sorted(chosen)})
    ind["ont"] = ind["sample"].isin(ont_ids)
    ind.to_csv(out / "individuals.tsv", sep="\t", index=False)
    r = runs[runs["sample"].isin(chosen)]
    r[["run", "sample", "read_pairs", "read_length", "sra"]].to_csv(out / "runs.tsv", sep="\t", index=False)
    have = M.loc[sorted(chosen)].sum()
    print(f"individuals: {len(chosen)} ({len(ont_ids)} with ONT, {len(extra)} added); runs: {len(r)}")
    print(f"testable genes (>= 10 informative in the whole cohort): {M.shape[1]}; with >= {args.min_informative} "
          f"informative in the chosen set: {int((have >= args.min_informative).sum())}")
    print(have.sort_values().to_string())


if __name__ == "__main__":
    main()
