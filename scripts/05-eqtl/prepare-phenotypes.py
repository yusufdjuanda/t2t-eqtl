"""Step 5: expression matrices -> tensorQTL phenotype BED for one cohort, reference and quantification method.

    pixi run python scripts/05-eqtl/prepare-phenotypes.py --cohort mage --ref chm13v2.0 [--method star|salmon]

Input: results/04-expression/<cohort>/matrices/<ref>/{star-counts,salmon-counts,salmon-tpm}.tsv.gz (merge-expression.py)
and data/star-references/<ref>/genes.tsv. GTEx-style processing:
  - genes on chr1-22 and chrX (T2T NC_ names converted to chrN with config/chrom-names-chm13v2.0.tsv, so the BED
    matches the genotype files); rRNA and chrM genes dropped
  - expressed: TPM >= 0.1 and >= 6 reads in >= 20% of individuals
  - TMM-normalised CPM (qtl.norm.edger_cpm), then inverse-normal transform per gene
  - one row per gene at its TSS (gene start on +, gene end on -); genes annotated at several loci keep their first
Writes results/05-eqtl/<cohort>/<ref>/<method>/phenotypes.bed.gz and genes-tested.tsv (gene_id, gene_ncbi, symbol).
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import qtl.norm

ROOT = Path(__file__).resolve().parents[2]
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--method", default="star", choices=["star", "salmon"])
    ap.add_argument("--tpm", type=float, default=0.1)
    ap.add_argument("--reads", type=int, default=6)
    ap.add_argument("--frac", type=float, default=0.2)
    args = ap.parse_args()

    mdir = ROOT / f"results/04-expression/{args.cohort}/matrices/{args.ref}"
    counts = pd.read_csv(mdir / f"{args.method}-counts.tsv.gz", sep="\t", index_col=0)
    tpm = pd.read_csv(mdir / "salmon-tpm.tsv.gz", sep="\t", index_col=0).reindex(index=counts.index, columns=counts.columns)
    genes = pd.read_csv(ROOT / f"data/star-references/{args.ref}/genes.tsv", sep="\t", dtype={"gene_ncbi": str})
    names = pd.read_csv(ROOT / "config/chrom-names-chm13v2.0.tsv", sep="\t")
    genes["chr"] = genes.chrom.map(dict(zip(names.refseq, names.ucsc))).fillna(genes.chrom)
    genes = genes[genes.chr.isin(CHROMS) & (genes.biotype != "rRNA")].drop_duplicates("gene_id")
    genes = genes.set_index("gene_id")

    counts = counts.loc[counts.index.intersection(genes.index)]
    tpm = tpm.loc[counts.index]
    keep = ((tpm >= args.tpm).mean(axis=1) >= args.frac) & ((counts >= args.reads).mean(axis=1) >= args.frac)
    print(f"{args.cohort} {args.ref} {args.method}: {counts.shape[1]} individuals, {len(counts)} genes on chr1-22/X, "
          f"{int(keep.sum())} expressed")
    # TMM over all genes (as GTEx), then keep the expressed ones
    cpm = qtl.norm.edger_cpm(counts, normalized_lib_sizes=True)
    norm = qtl.norm.inverse_normal_transform(cpm.loc[keep])

    g = genes.loc[norm.index]
    tss = np.where(g.strand == "+", g.start, g.end)
    bed = pd.DataFrame({"chr": g.chr, "start": tss - 1, "end": tss, "phenotype_id": g.index})
    bed = pd.concat([bed.reset_index(drop=True), norm.reset_index(drop=True)], axis=1)
    bed["_c"] = bed.chr.map({c: i for i, c in enumerate(CHROMS)})
    bed = bed.sort_values(["_c", "start"]).drop(columns="_c")

    out = ROOT / f"results/05-eqtl/{args.cohort}/{args.ref}/{args.method}"
    out.mkdir(parents=True, exist_ok=True)
    bed.rename(columns={"chr": "#chr"}).to_csv(out / "phenotypes.bed.gz", sep="\t", index=False, float_format="%.5g")
    g.reset_index()[["gene_id", "gene_ncbi", "symbol", "chr", "start", "end", "strand"]].to_csv(
        out / "genes-tested.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
