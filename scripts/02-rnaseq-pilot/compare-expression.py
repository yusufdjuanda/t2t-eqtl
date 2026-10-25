"""Step 2 (go/no-go 1): compare pilot gene expression between T2T-CHM13v2.0 and GRCh38.

    pixi run python scripts/02-rnaseq-pilot/compare-expression.py [--sheet config/samples/pilot-5.tsv]

Reads results/02-rnaseq-pilot/<ref>/<run>/ (STAR ReadsPerGene, salmon quant.sf, Log.final.out, idxstats) and
data/star-references/<ref>/genes.tsv, tx2gene.tsv. Genes are matched across references by NCBI GeneID.
Writes results/02-rnaseq-pilot/comparison/:
  mapping-stats.tsv   per run x reference: unique / multi / unmapped %, chrM and chrEBV share of mapped records
  genes.tsv           per GeneID: loci on each reference, T2T segmental-duplication overlap, counts and CPM per
                      sample x reference x method, mean log2(T2T/GRCh38) and sign consistency per method
  candidates.tsv      the go/no-go candidate families (chr21 false duplication, 1q21, 5q13, HYDIN)
  summary.tsv         expressed and changed genes per method, and the SD share among changed vs all expressed
  scatter-<method>.png
CPM uses genes on both references minus rRNA and chrM genes as the library size.
A gene is "changed" for a method when mean CPM >= 1 on either reference, |mean log2 ratio| >= 1, and all
samples change in the same direction.
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
REFS = ["chm13v2.0", "grch38"]
METHODS = ["star", "salmon"]
PSEUDO = 0.5  # CPM pseudocount for log ratios
CANDIDATES = {  # family: symbols; LOC… are the GRCh38-only false copies on chr21p
    "chr21 CBS": ["CBS", "LOC102724560"],
    "chr21 U2AF1": ["U2AF1", "LOC102724594"],
    "chr21 CRYAA": ["CRYAA", "LOC102724652"],
    "chr21 SIK1": ["SIK1", "LOC102724428"],
    "chr21 ICOSLG": ["ICOSLG", "LOC102723996"],
    "chr21 DNMT3L": ["DNMT3L", "LOC124905051"],
    "chr21 GATD3": ["GATD3", "LOC102724023"],
    "chr21 PWP2": ["PWP2", "LOC102724159"],
    "chr21 TRAPPC10": ["TRAPPC10", "LOC102724200"],
    "chr21 KCNE1": ["KCNE1"],
    "1q21 NOTCH2NL": ["NOTCH2", "NOTCH2NLA", "NOTCH2NLB", "NOTCH2NLC", "NOTCH2NLR"],
    "1q21 NBPF": ["NBPF*"],
    "5q13 SMN": ["SMN1", "SMN2", "NAIP", "SERF1A", "SERF1B", "GTF2H2", "GTF2H2C"],
    "16q22 HYDIN": ["HYDIN", "HYDIN2"],
}


def gene_table(ref):
    g = pd.read_csv(ROOT / f"data/star-references/{ref}/genes.tsv", sep="\t", dtype={"gene_ncbi": str})
    return g


def star_counts(ref, run, column="unstranded"):
    t = pd.read_csv(ROOT / f"results/02-rnaseq-pilot/{ref}/{run}/ReadsPerGene.out.tab", sep="\t", header=None,
                    names=["gene_id", "unstranded", "fwd", "rev"])
    return t[~t.gene_id.str.startswith("N_")].set_index("gene_id")[column]


def salmon_counts(ref, run, tx2gene):
    q = pd.read_csv(ROOT / f"results/02-rnaseq-pilot/{ref}/{run}/salmon/quant.sf", sep="\t")
    q["gene_id"] = q.Name.map(tx2gene)
    return q.groupby("gene_id").NumReads.sum()


def log_final(ref, run):
    out = {}
    for line in open(ROOT / f"results/02-rnaseq-pilot/{ref}/{run}/Log.final.out"):
        if "|" in line:
            k, v = (x.strip() for x in line.split("|"))
            out[k] = v
    idx = pd.read_csv(ROOT / f"results/02-rnaseq-pilot/{ref}/{run}/idxstats.tsv", sep="\t", header=None,
                      names=["contig", "length", "mapped", "unmapped"]).set_index("contig")
    mapped = idx.mapped.sum()
    pct = lambda k: float(out[k].rstrip("%"))
    return {
        "input_pairs": int(out["Number of input reads"]),
        "unique_pct": pct("Uniquely mapped reads %"),
        "multi_pct": pct("% of reads mapped to multiple loci"),
        "too_many_loci_pct": pct("% of reads mapped to too many loci"),
        "unmapped_too_short_pct": pct("% of reads unmapped: too short"),
        "unmapped_other_pct": pct("% of reads unmapped: other"),  # on T2T mostly rRNA (seeds hit too many rDNA copies)
        "chrM_pct_of_mapped_records": 100 * idx.mapped.get("chrM", 0) / mapped,
        "chrEBV_pct_of_mapped_records": 100 * idx.mapped.get("chrEBV", 0) / mapped,
    }


def sd_overlap(genes):
    """Fraction of each T2T gene covered by segmental duplications (chm13v2.0_SD.bed, UCSC names)."""
    names = pd.read_csv(ROOT / "config/chrom-names-chm13v2.0.tsv", sep="\t")
    to_refseq = dict(zip(names.ucsc, names.refseq))
    sd = pd.read_csv(ROOT / "data/references/t2t-annotations/chm13v2.0_SD.bed", sep="\t", header=None, usecols=[0, 1, 2],
                     names=["chrom", "start", "end"])
    sd["chrom"] = sd.chrom.map(to_refseq)
    frac = pd.Series(0.0, index=genes.index)
    for chrom, iv in sd.groupby("chrom"):
        iv = iv.sort_values("start")
        merged = []  # merge overlapping SD intervals (pairs are listed from both sides)
        for s, e in zip(iv.start, iv.end):
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        m = np.array(merged)
        cum = np.concatenate([[0], np.cumsum(m[:, 1] - m[:, 0])])

        def covered(x):  # bases of [0, x) covered by merged intervals
            i = np.searchsorted(m[:, 0], x, side="right")
            prev = cum[i - 1] if i > 0 else 0
            return prev + np.clip(x - m[i - 1, 0], 0, m[i - 1, 1] - m[i - 1, 0]) if i > 0 else 0

        g = genes[genes.chrom == chrom]
        for idx, s, e in zip(g.index, g.start - 1, g.end):
            frac[idx] = (covered(e) - covered(s)) / (e - s)
    return frac


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sheet", default=ROOT / "config/samples/pilot-5.tsv")
    ap.add_argument("--out", default="comparison", help="output folder name under results/02-rnaseq-pilot/")
    ap.add_argument("--star-column", default="unstranded", choices=["unstranded", "fwd", "rev"],
                    help="ReadsPerGene column: unstranded (Geuvadis) or the strand matching a stranded library")
    args = ap.parse_args()
    sheet = pd.read_csv(args.sheet, sep="\t")
    runs = list(sheet.run)
    out = ROOT / "results/02-rnaseq-pilot" / args.out
    out.mkdir(parents=True, exist_ok=True)

    stats = pd.DataFrame([{"run": r, "sample": s, "reference": ref, **log_final(ref, r)}
                          for r, s in zip(sheet.run, sheet["sample"]) for ref in REFS])
    stats.to_csv(out / "mapping-stats.tsv", sep="\t", index=False, float_format="%.3f")

    genes = {ref: gene_table(ref) for ref in REFS}
    genes["chm13v2.0"]["sd_frac"] = sd_overlap(genes["chm13v2.0"])
    table = None
    for ref in REFS:
        g = genes[ref]
        tx2gene = pd.read_csv(ROOT / f"data/star-references/{ref}/tx2gene.tsv", sep="\t")
        tx2gene = dict(zip(tx2gene.transcript_id, tx2gene.gene_id))
        cols = {}
        for run in runs:
            cols[f"star_{ref}_{run}"] = star_counts(ref, run, args.star_column)
            cols[f"salmon_{ref}_{run}"] = salmon_counts(ref, run, tx2gene)
        c = pd.DataFrame(cols).reindex(g.gene_id).fillna(0)
        c["gene_ncbi"] = g.gene_ncbi.values
        agg = c.groupby("gene_ncbi").sum()  # a GeneID annotated at several loci is summed
        locus = g.assign(locus=g.chrom + ":" + g.start.astype(str) + "-" + g.end.astype(str)).groupby("gene_ncbi").agg(
            symbol=("symbol", "first"), biotype=("biotype", "first"), n_loci=("locus", "size"),
            loci=("locus", lambda x: ";".join(x)), **({"sd_frac": ("sd_frac", "max")} if ref == "chm13v2.0" else {}))
        locus = locus.rename(columns={k: f"{k}_{ref}" for k in ["symbol", "biotype", "n_loci", "loci"]})
        part = locus.join(agg)
        table = part if table is None else table.join(part, how="outer")

    table["symbol"] = table["symbol_chm13v2.0"].fillna(table["symbol_grch38"])
    table["biotype"] = table["biotype_chm13v2.0"].fillna(table["biotype_grch38"])
    for ref in REFS:
        table[f"n_loci_{ref}"] = table[f"n_loci_{ref}"].fillna(0).astype(int)
    count_cols = [c for c in table.columns if c.startswith(("star_", "salmon_"))]
    table[count_cols] = table[count_cols].fillna(0)

    # Library size for CPM: genes with a GeneID on both references, minus rRNA and chrM genes. T2T has the full rDNA
    # arrays (STAR leaves those reads unmapped; only 14 45S-type genes are annotated), while salmon gives GRCh38's five
    # rDNA copies ~10,000 CPM; counting either in the denominator would shift every other gene between references.
    both = (table["n_loci_chm13v2.0"] > 0) & (table["n_loci_grch38"] > 0)
    norm_genes = both & (table.biotype != "rRNA") & ~table["loci_grch38"].fillna("").str.startswith("chrM:")
    summary = []
    for method in METHODS:
        for ref in REFS:
            cs = [f"{method}_{ref}_{r}" for r in runs]
            cpm = table[cs] / table.loc[norm_genes, cs].sum() * 1e6
            table[[f"cpm_{c}" for c in cs]] = cpm.values
            table[f"mean_cpm_{method}_{ref}"] = cpm.mean(axis=1)
        lr = np.column_stack([np.log2((table[f"cpm_{method}_chm13v2.0_{r}"] + PSEUDO) /
                                      (table[f"cpm_{method}_grch38_{r}"] + PSEUDO)) for r in runs])
        table[f"log2_t2t_vs_grch38_{method}"] = lr.mean(axis=1)
        table[f"same_sign_{method}"] = (np.all(lr > 0, axis=1) | np.all(lr < 0, axis=1))
        expressed = table[[f"mean_cpm_{method}_{ref}" for ref in REFS]].max(axis=1) >= 1
        changed = expressed & table[f"same_sign_{method}"] & (table[f"log2_t2t_vs_grch38_{method}"].abs() >= 1)
        table[f"changed_{method}"] = changed
        in_sd = table.sd_frac.fillna(0) >= 0.5
        summary.append({
            "method": method, "expressed": int(expressed.sum()),
            "expressed_both_refs": int((expressed & both).sum()),
            "expressed_t2t_only_gene": int((expressed & (table["n_loci_grch38"] == 0)).sum()),
            "expressed_grch38_only_gene": int((expressed & (table["n_loci_chm13v2.0"] == 0)).sum()),
            "changed_up_on_t2t": int((changed & both & (table[f"log2_t2t_vs_grch38_{method}"] > 0)).sum()),
            "changed_down_on_t2t": int((changed & both & (table[f"log2_t2t_vs_grch38_{method}"] < 0)).sum()),
            "sd_pct_of_expressed": round(100 * (in_sd & expressed & both).sum() / (expressed & both).sum(), 2),
            "sd_pct_of_changed": round(100 * (in_sd & changed & both).sum() / max((changed & both).sum(), 1), 2),
            "median_abs_log2_expressed": round(table.loc[expressed & both, f"log2_t2t_vs_grch38_{method}"].abs().median(), 4),
        })

        fig, ax = plt.subplots(figsize=(5.5, 5.5))
        x = np.log10(table.loc[both, f"mean_cpm_{method}_grch38"] + PSEUDO)
        y = np.log10(table.loc[both, f"mean_cpm_{method}_chm13v2.0"] + PSEUDO)
        sd_mask = in_sd[both]
        ax.scatter(x[~sd_mask], y[~sd_mask], s=2, alpha=0.3, color="0.5", label="other genes", rasterized=True)
        ax.scatter(x[sd_mask], y[sd_mask], s=4, alpha=0.6, color="tab:red", label="≥50% in T2T SDs", rasterized=True)
        lim = [np.log10(PSEUDO), max(x.max(), y.max())]
        ax.plot(lim, lim, color="k", lw=0.5)
        ax.set_xlabel("GRCh38 mean CPM (log10, +0.5)")
        ax.set_ylabel("T2T-CHM13 mean CPM (log10, +0.5)")
        ax.set_title(f"{method}: {len(runs)} Geuvadis samples, genes on both references")
        ax.legend(loc="upper left", markerscale=3, frameon=False)
        fig.tight_layout()
        fig.savefig(out / f"scatter-{method}.png", dpi=150)
        plt.close(fig)

    table.index.name = "gene_ncbi"
    table.to_csv(out / "genes.tsv", sep="\t", float_format="%.4g")
    pd.DataFrame(summary).to_csv(out / "summary.tsv", sep="\t", index=False)

    rows = []
    for family, symbols in CANDIDATES.items():
        for sym in symbols:
            hit = table[table.symbol.str.startswith(sym[:-1])] if sym.endswith("*") else table[table.symbol == sym]
            for gid, r in hit.iterrows():
                rows.append({"family": family, "symbol": r.symbol, "gene_ncbi": gid,
                             "n_loci_chm13v2.0": r["n_loci_chm13v2.0"], "n_loci_grch38": r["n_loci_grch38"],
                             **{f"mean_cpm_{m}_{ref}": r[f"mean_cpm_{m}_{ref}"] for m in METHODS for ref in REFS},
                             **{f"log2_t2t_vs_grch38_{m}": r[f"log2_t2t_vs_grch38_{m}"] for m in METHODS},
                             **{f"same_sign_{m}": r[f"same_sign_{m}"] for m in METHODS}})
    cand = pd.DataFrame(rows)
    cand.to_csv(out / "candidates.tsv", sep="\t", index=False, float_format="%.3f")

    pd.set_option("display.width", 250)
    print(stats.to_string(index=False, float_format="%.2f"))
    print(pd.DataFrame(summary).to_string(index=False))
    print(cand[["family", "symbol", "n_loci_chm13v2.0", "n_loci_grch38", "mean_cpm_star_chm13v2.0", "mean_cpm_star_grch38",
                "log2_t2t_vs_grch38_star", "mean_cpm_salmon_chm13v2.0", "mean_cpm_salmon_grch38",
                "log2_t2t_vs_grch38_salmon"]].to_string(index=False, float_format="%.2f"))


if __name__ == "__main__":
    main()
