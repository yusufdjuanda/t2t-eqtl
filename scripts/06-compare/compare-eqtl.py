"""Step 6: compare cis-eQTL results between references (T2T-CHM13, GRCh38, FixItFelix-patched GRCh38) per gene.

    pixi run python scripts/06-compare/compare-eqtl.py --cohort mage420 [--method star] [--variants matched] [--fdr 0.05]

Inputs, per reference: results/05-eqtl/<cohort>/<ref>/<method>/{genes-tested.tsv, <variants>/cis-eqtl.tsv.gz,
<variants>/nominal/*.parquet}; results/06-compare/{gene-classes.tsv, annotation-overlaps.tsv, variant-matching/}.
Genes are matched by NCBI GeneID; lead variants are translated between T2T and GRCh38 IDs with the variant-matching
tables (both GRCh38 references share GRCh38 IDs).
Per gene and reference pair (A, B) the table records eGene status (q <= fdr), and for an eGene on A whether its lead
variant also associates on B (nominal p < 0.05 for the same gene and the matched variant): "replicates on B".
A gene counts as a reference-specific eGene for A vs B ("A_only") only when it is an eGene on A, not on B, and its
lead variant does not replicate on B (or the gene / variant is untestable on B) - so threshold effects are not counted.
Writes results/06-compare/eqtl-comparison/<cohort>/<method>-<variants>/genes.tsv and summary.tsv, and prints the
H0 summary: reference-specific eGenes by region class, without artefact flags (annotation overlap differs, MHC, PAR).
"""
import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
REFS = ["chm13v2.0", "grch38", "grch38-ffx"]
SIDE = {"chm13v2.0": "t2t", "grch38": "grch38", "grch38-ffx": "grch38"}


def load_ref(cohort, ref, method, variants):
    d = ROOT / f"results/05-eqtl/{cohort}/{ref}/{method}"
    genes = pd.read_csv(d / "genes-tested.tsv", sep="\t", dtype={"gene_ncbi": str})
    cis = pd.read_csv(d / variants / "cis-eqtl.tsv.gz", sep="\t", index_col=0)
    cis = cis.join(genes.set_index("gene_id")[["gene_ncbi"]])
    return cis


def nominal_lookup(cohort, ref, method, variants, pairs):
    """p-values for (phenotype_id, variant_id) pairs from the nominal parquet files of one reference."""
    if pairs.empty:
        return pd.Series(dtype=float)
    want = set(map(tuple, pairs[["phenotype_id", "variant_id"]].values))
    genes = set(pairs.phenotype_id)
    out = {}
    for f in glob.glob(str(ROOT / f"results/05-eqtl/{cohort}/{ref}/{method}/{variants}/nominal/*.parquet")):
        x = pd.read_parquet(f, columns=["phenotype_id", "variant_id", "pval_nominal"])
        x = x[x.phenotype_id.isin(genes)]
        for g, v, p in zip(x.phenotype_id, x.variant_id, x.pval_nominal):
            if (g, v) in want:
                out[(g, v)] = p
    return pd.Series(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--method", default="star")
    ap.add_argument("--variants", default="matched")
    ap.add_argument("--fdr", type=float, default=0.05)
    args = ap.parse_args()

    # variant ID translation T2T <-> GRCh38 (matched pairs only)
    m = pd.concat(pd.read_csv(f, sep="\t", usecols=["side", "category", "t2t_id", "grch38_id"], low_memory=False)
                  for f in sorted(glob.glob(str(ROOT / "results/06-compare/variant-matching/chr*.tsv.gz"))))
    m = m[(m.side == "t2t") & (m.category == "matched")]
    t2g = dict(zip(m.t2t_id, m.grch38_id)); g2t = dict(zip(m.grch38_id, m.t2t_id))

    def translate(v, a, b):
        if SIDE[a] == SIDE[b]:
            return v
        return t2g.get(v) if SIDE[a] == "t2t" else g2t.get(v)

    res = {r: load_ref(args.cohort, r, args.method, args.variants) for r in REFS}
    # one row per GeneID; a GeneID tested at several loci keeps its best gene
    per = {}
    for r, c in res.items():
        c = c.dropna(subset=["gene_ncbi"]).sort_values("pval_beta").drop_duplicates("gene_ncbi")
        per[r] = c.reset_index().set_index("gene_ncbi")[["phenotype_id", "variant_id", "slope", "pval_beta", "qval"]]
    t = pd.concat({r: per[r] for r in REFS}, axis=1)
    t.columns = [f"{c}_{r}" for r, c in t.columns]
    for r in REFS:
        t[f"tested_{r}"] = t[f"qval_{r}"].notna()
        t[f"egene_{r}"] = t[f"qval_{r}"] <= args.fdr

    # replication of each eGene's lead variant on the other references
    for a in REFS:
        for b in REFS:
            if a == b:
                continue
            e = t[t[f"egene_{a}"] & t[f"tested_{b}"]]
            pairs = pd.DataFrame({"gene": e.index, "phenotype_id": e[f"phenotype_id_{b}"].values,
                                  "variant_id": [translate(v, a, b) for v in e[f"variant_id_{a}"]]}).dropna()
            p = nominal_lookup(args.cohort, b, args.method, args.variants, pairs)
            pv = pd.Series([p.get((g, v), np.nan) for g, v in zip(pairs.phenotype_id, pairs.variant_id)], index=pairs.gene)
            t[f"rep_{a}_on_{b}"] = pv.reindex(t.index) < 0.05
            t[f"only_{a}_vs_{b}"] = t[f"egene_{a}"] & ~t[f"egene_{b}"] & ~t[f"rep_{a}_on_{b}"].fillna(False).astype(bool)

    cls = pd.read_csv(ROOT / "results/06-compare/gene-classes.tsv", sep="\t", dtype={"gene_ncbi": str}).set_index("gene_ncbi")
    ovl = pd.read_csv(ROOT / "results/06-compare/annotation-overlaps.tsv", sep="\t", dtype={"gene_ncbi": str}).set_index("gene_ncbi")
    t = t.join(cls[["symbol", "biotype", "annotated_on", "region_class", "flag_mhc", "flag_par"]]).join(
        ovl[["overlap_differs"]])
    t["overlap_differs"] = t.overlap_differs.astype(str).eq("True")
    t["artefact_flag"] = t.overlap_differs | t.flag_mhc.fillna(False).astype(bool) | t.flag_par.fillna(False).astype(bool)

    out = ROOT / f"results/06-compare/eqtl-comparison/{args.cohort}/{args.method}-{args.variants}"
    out.mkdir(parents=True, exist_ok=True)
    t.index.name = "gene_ncbi"
    t.to_csv(out / "genes.tsv", sep="\t", float_format="%.4g")

    rows = []
    for r in REFS:
        rows.append({"comparison": f"eGenes on {r}", "all": int(t[f"egene_{r}"].sum())})
    for a, b in [("chm13v2.0", "grch38"), ("grch38", "chm13v2.0"), ("grch38-ffx", "grch38"), ("grch38", "grch38-ffx"),
                 ("chm13v2.0", "grch38-ffx"), ("grch38-ffx", "chm13v2.0")]:
        x = t[t[f"only_{a}_vs_{b}"]]
        clean = x[~x.artefact_flag]
        row = {"comparison": f"{a} only vs {b}", "all": len(x), "no_artefact_flag": len(clean)}
        row.update({f"clean_{k}": int(v) for k, v in clean.region_class.value_counts().items()})
        rows.append(row)
    s = pd.DataFrame(rows).fillna(0)
    s.to_csv(out / "summary.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(s.to_string(index=False))


if __name__ == "__main__":
    main()
