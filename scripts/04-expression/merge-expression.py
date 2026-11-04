"""Stage 04 (expression): merge per-run quantifications into gene x individual matrices, one set per reference.

    pixi run python scripts/04-expression/merge-expression.py --cohort mage [--refs chm13v2.0,grch38,grch38-ffx] [--jobs 3]

Reads results/04-expression/<cohort>/<run>/ for every run in config/samplesheets/<cohort>.csv that has complete.txt.
Writes results/04-expression/<cohort>/matrices/:
  runs.tsv                         every samplesheet run: included or not, and its mapping summary per reference
  <ref>/star-counts.tsv.gz         STAR GeneCounts, unstranded column (uniquely mapped read pairs), summed over an
                                   individual's runs
  <ref>/salmon-counts.tsv.gz       salmon NumReads summed to genes (tx2gene), summed over runs
  <ref>/salmon-tpm.tsv.gz          salmon TPM summed to genes, averaged over an individual's runs
Rows: gene_id of that reference's genes.tsv (match across references by gene_ncbi there). Columns: individuals.
Both MAGE and Geuvadis are unstranded (checked in the pilots), hence the unstranded STAR column.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def star(path):
    t = pd.read_csv(path, sep="\t", header=None, names=["gene_id", "unstranded", "fwd", "rev"])
    return t[~t.gene_id.str.startswith("N_")].set_index("gene_id")["unstranded"]


def salmon(path, tx2gene):
    q = pd.read_csv(path, sep="\t", usecols=["Name", "TPM", "NumReads"])
    q["gene_id"] = q.Name.map(tx2gene)
    g = q.groupby("gene_id")[["NumReads", "TPM"]].sum()
    return g.NumReads, g.TPM


def merge_ref(ref, cohort, runs):
    res = ROOT / "results/04-expression" / cohort
    genes = pd.read_csv(ROOT / f"data/star-references/{ref}/genes.tsv", sep="\t").gene_id
    tx = pd.read_csv(ROOT / f"data/star-references/{ref}/tx2gene.tsv", sep="\t")
    tx2gene = dict(zip(tx.transcript_id, tx.gene_id))
    sc, nr, tpm = {}, {}, {}
    for run in runs.run:
        sc[run] = star(res / run / ref / "ReadsPerGene.out.tab")
        nr[run], tpm[run] = salmon(res / run / ref / "salmon/quant.sf", tx2gene)
    by = runs.set_index("run")["sample"]
    out = res / "matrices" / ref
    out.mkdir(parents=True, exist_ok=True)
    for name, d, how in [("star-counts", sc, "sum"), ("salmon-counts", nr, "sum"), ("salmon-tpm", tpm, "mean")]:
        m = pd.DataFrame(d).reindex(genes).fillna(0)
        m = m.T.groupby(by).agg(how).T  # runs -> individuals
        m.index.name = "gene_id"
        m.round(3).to_csv(out / f"{name}.tsv.gz", sep="\t", compression="gzip")
    return ref, m.shape


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--refs", default="chm13v2.0,grch38,grch38-ffx")
    ap.add_argument("--jobs", type=int, default=3)
    args = ap.parse_args()
    refs = args.refs.split(",")
    res = ROOT / "results/04-expression" / args.cohort
    sheet = pd.read_csv(ROOT / f"config/samplesheets/{args.cohort}.csv")
    sheet["included"] = [(res / r / "complete.txt").exists() for r in sheet.run]
    summ = [pd.read_csv(res / r / "mapping-summary.tsv", sep="\t") for r in sheet.run[sheet.included]]
    if summ:
        wide = pd.concat(summ).pivot(index="run", columns="reference",
                                     values=["unique_pct", "multi_pct", "unmapped_other_pct"])
        wide.columns = [f"{v}_{r}" for v, r in wide.columns]
        sheet = sheet.merge(wide, left_on="run", right_index=True, how="left")
    (res / "matrices").mkdir(parents=True, exist_ok=True)
    sheet.drop(columns=["sra"]).to_csv(res / "matrices/runs.tsv", sep="\t", index=False)
    runs = sheet[sheet.included][["run", "sample"]]
    print(f"{args.cohort}: {len(runs)} of {len(sheet)} runs complete, {runs['sample'].nunique()} individuals")
    if runs.empty:
        return
    with ProcessPoolExecutor(args.jobs) as ex:
        for ref, shape in ex.map(merge_ref, refs, [args.cohort] * len(refs), [runs] * len(refs)):
            print(f"  {ref}: {shape[0]} genes x {shape[1]} individuals")


if __name__ == "__main__":
    main()
