"""Step 6: per-gene paralog copy numbers of the 1KGP individuals from ctyper (Ma & Chaisson 2025), for use as covariates.

    pixi run python scripts/06-compare/ctyper-copy-number.py --table <PangenomeAlleles_typefix.tsv_1KG_summary.txt>

Input: ctyper's 1KGP genotype table (data/external/ctyper-ma2025/PangenomeAlleles_1kg.tar.gz, extracted to scratch):
one row per pangenome allele type (84,346), columns 1-10 annotation, then the copy number of that allele type in each
of 3,202 individuals. Column 1 is <family>_group<k>_<i>_<j>; column 3 lists the genes the allele overlaps
("A,B:C,D(part1)&..."); column 4 is the main gene, "<symbol>:<count>".
Per gene symbol of our three annotations (data/star-references/*/genes.tsv):
  gene     sum of the allele types whose main gene is the symbol (copies called as that gene)
  family   sum over all allele types of every ctyper group that mentions the symbol in column 3 or 4 (all copies
           of the paralog family whose reads can be confused with it)
Symbols are matched as written; ctyper uses GENCODE names on GRCh38, so some RefSeq LOC genes are absent.
Writes results/06-compare/ctyper/copy-number.tsv.gz (symbol, kind, then one column per individual of MAGE and
Geuvadis) and prints coverage.
"""
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    args = ap.parse_args()
    out = ROOT / "results/06-compare/ctyper"
    out.mkdir(parents=True, exist_ok=True)
    keep = set()
    for c in ("mage", "geuvadis"):
        keep |= set(pd.read_csv(ROOT / f"results/04-expression/{c}/matrices/runs.tsv", sep="\t")["sample"])
    symbols = set()
    for r in ("chm13v2.0", "grch38", "grch38-ffx"):
        symbols |= set(pd.read_csv(ROOT / f"data/star-references/{r}/genes.tsv", sep="\t").symbol.dropna())

    header = pd.read_csv(args.table, sep="\t", nrows=0).columns
    samples = [h.split(".")[0] for h in header[1:]]
    ncol = len(pd.read_csv(args.table, sep="\t", nrows=1, header=None, skiprows=1).columns)
    nann = ncol - len(samples)
    want = [i for i, s in enumerate(samples) if s in keep]
    t = pd.read_csv(args.table, sep="\t", header=None, skiprows=1, usecols=list(range(nann)) + [nann + i for i in want])
    ann = t.iloc[:, :nann]
    cn = t.iloc[:, nann:].to_numpy(dtype=float)
    ids = [samples[i] for i in want]
    print(f"{len(t)} allele types, {nann} annotation columns, {len(ids)} of {len(keep)} cohort individuals present")

    group = ann[0].str.extract(r"^(.*_group\d+)_")[0].fillna(ann[0])
    main = ann[3].astype(str).str.split(":").str[0]
    tokens = [set(x for x in re.split(r"[,:&;]", re.sub(r"\(part\d+\)", "", f"{g}:{m}")) if x)
              for g, m in zip(ann[2].astype(str), main)]
    gidx = {g: np.flatnonzero(group.values == g) for g in pd.unique(group)}
    sym_groups = {}
    for tk, g in zip(tokens, group):
        for s in tk & symbols:
            sym_groups.setdefault(s, set()).add(g)
    rows, vals = [], []
    for s, idx in pd.Series(np.arange(len(t))).groupby(main.values):
        if s in symbols:
            rows.append((s, "gene")); vals.append(cn[idx.values].sum(0))
    for s, gs in sym_groups.items():
        idx = np.concatenate([gidx[g] for g in gs])
        rows.append((s, "family")); vals.append(cn[idx].sum(0))
    res = pd.DataFrame(np.vstack(vals), columns=ids)
    res.insert(0, "kind", [k for _, k in rows]); res.insert(0, "symbol", [s for s, _ in rows])
    res = res.sort_values(["symbol", "kind"])
    res.to_csv(out / "copy-number.tsv.gz", sep="\t", index=False, float_format="%.3g")
    var = res.iloc[:, 2:].var(axis=1) > 0
    print(f"symbols with a gene copy number: {int((res.kind == 'gene').sum())}, with a family copy number: "
          f"{int((res.kind == 'family').sum())}; variable across individuals: {int(var.sum())} of {len(res)} rows")


if __name__ == "__main__":
    main()
