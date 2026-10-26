"""Presentation: export small tables for the R figures (make-figures.R) and the workbook (build-workbook.py).

    pixi run python scripts/07-presentation/prepare-figure-data.py

Reads finished results only (MAGE; Geuvadis appears as a cohort, its expression is still running). Writes
results/07-presentation/data/*.tsv:
  cohorts.tsv           one row per run of both samplesheets (cohort, population, superpopulation, sex, read pairs)
  mapping-qc.tsv        MAGE runs.tsv: uniquely / multi-mapped / unmapped % per run and reference
  variant-matching.tsv  results/06-compare/variant-matching/summary.tsv
  egenes.tsv            genes tested and eGenes per reference and quantification method
  effects.tsv           same-variant comparison T2T vs patched GRCh38 (STAR), all genes, with robust flag
  positive-control.tsv  q-values on the three references for the FixItFelix genes
  robust-genes.tsv      the robust set joined with its read destination
  paralog-deltas.tsv    per-individual CPM change (T2T - patched GRCh38) for example gene pairs
  eqtl-examples.tsv     normalised expression and lead-variant genotype on both references for example genes
"""
import glob
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/07-presentation/data"
CMP = ROOT / "results/06-compare/eqtl-comparison/mage"
REFS = ["chm13v2.0", "grch38", "grch38-ffx"]
PAIRS = [("CCZ1", "CCZ1B"), ("ASAH2", "ASAH2B"), ("RPSA", "RPSA2"), ("PDPR", "LOC128966601")]
EQTL_EXAMPLES = ["CCZ1", "SMN2", "PWP2"]


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    # cohorts: population codes look like "MXL (AMR)" in both samplesheets
    sheets = [pd.read_csv(ROOT / f"config/samplesheets/{c}.csv") for c in ("mage", "geuvadis")]
    c = pd.concat(sheets)[["cohort", "sample", "run", "ancestry", "sex", "read_pairs", "read_length", "lab"]]
    # MAGE writes "MXL (AMR)"; Geuvadis (E-GEUV-1) writes population names
    geuv = {"Tuscan": ("TSI", "EUR"), "Utah": ("CEU", "EUR"), "Finnish": ("FIN", "EUR"), "British": ("GBR", "EUR"),
            "Yoruba": ("YRI", "AFR")}
    pop = c.ancestry.str.extract(r"^(?P<population>[^ (]+)(?: \((?P<superpop>[A-Z]+)\))?")
    named = pop.population.isin(geuv)
    pop.loc[named, "superpop"] = pop.loc[named, "population"].map(lambda p: geuv[p][1])
    pop.loc[named, "population"] = pop.loc[named, "population"].map(lambda p: geuv[p][0])
    c = pd.concat([c.reset_index(drop=True), pop.reset_index(drop=True)], axis=1)
    c.to_csv(OUT / "cohorts.tsv", sep="\t", index=False)

    pd.read_csv(ROOT / "results/03-expression/mage/matrices/runs.tsv", sep="\t").to_csv(OUT / "mapping-qc.tsv", sep="\t", index=False)
    pd.read_csv(ROOT / "results/06-compare/variant-matching/summary.tsv", sep="\t").to_csv(OUT / "variant-matching.tsv", sep="\t", index=False)

    rows = []
    for m in ("star", "salmon"):
        for r in REFS:
            x = pd.read_csv(ROOT / f"results/05-eqtl/mage/{r}/{m}/matched/cis-eqtl.tsv.gz", sep="\t")
            rows.append({"reference": r, "method": m, "genes_tested": len(x), "egenes": int((x.qval <= 0.05).sum())})
    pd.DataFrame(rows).to_csv(OUT / "egenes.tsv", sep="\t", index=False)

    e = pd.read_csv(CMP / "star-matched/effects-chm13v2.0-vs-grch38-ffx.tsv", sep="\t", dtype={"gene_ncbi": str})
    robust = pd.read_csv(CMP / "robust-matched/changed-chm13v2.0-vs-grch38-ffx.tsv", sep="\t", dtype={"gene_ncbi": str})
    e["robust"] = e.gene_ncbi.isin(robust.loc[~robust.artefact_flag, "gene_ncbi"])
    e.to_csv(OUT / "effects.tsv", sep="\t", index=False)

    g = pd.read_csv(CMP / "star-matched/genes.tsv", sep="\t", dtype={"gene_ncbi": str}, low_memory=False)
    pc = g[g.symbol.isin(["PWP2", "GATD3", "CBS", "U2AF1", "ICOSLG"])]
    pc = pc.melt(id_vars=["symbol"], value_vars=[f"qval_{r}" for r in REFS], var_name="reference", value_name="qval")
    pc["reference"] = pc.reference.str.replace("qval_", "")
    pc.to_csv(OUT / "positive-control.tsv", sep="\t", index=False)

    rd = pd.read_csv(CMP / "robust-matched/read-destination-chm13v2.0-vs-grch38-ffx.tsv", sep="\t", dtype={"gene_ncbi": str})
    rg = robust[~robust.artefact_flag].merge(rd.drop(columns=["symbol", "region_class", "stronger_on"]), on="gene_ncbi", how="left")
    rg.to_csv(OUT / "robust-genes.tsv", sep="\t", index=False)

    # per-individual CPM deltas for example pairs (same construction as read-destination.py)
    cpm = {}
    for r in ("chm13v2.0", "grch38-ffx"):
        m = pd.read_csv(ROOT / f"results/03-expression/mage/matrices/{r}/star-counts.tsv.gz", sep="\t", index_col=0)
        gt = pd.read_csv(ROOT / f"data/star-references/{r}/genes.tsv", sep="\t")
        m = m.groupby(m.index.map(dict(zip(gt.gene_id, gt.symbol)))).sum()
        cpm[r] = m / m.sum() * 1e6
    ids = cpm["chm13v2.0"].columns.intersection(cpm["grch38-ffx"].columns)
    rows = []
    for a, b in PAIRS:
        for gene, role in ((a, "gene"), (b, "partner")):
            ca = cpm["chm13v2.0"].reindex([gene], columns=ids).fillna(0).iloc[0]
            cb = cpm["grch38-ffx"].reindex([gene], columns=ids).fillna(0).iloc[0]
            rows.append(pd.DataFrame({"pair": f"{a} / {b}", "role": role, "symbol": gene, "sample": ids,
                                      "cpm_t2t": ca.values, "cpm_ffx": cb.values, "delta": (ca - cb).values}))
    pd.concat(rows).to_csv(OUT / "paralog-deltas.tsv", sep="\t", index=False)

    # eQTL examples: normalised expression vs lead-variant genotype (T2T lead, its matched GRCh38 variant)
    from tensorqtl import pgen
    ga = pgen.PgenReader(str(ROOT / "results/05-eqtl/mage/genotypes/t2t.matched"))
    gb = pgen.PgenReader(str(ROOT / "results/05-eqtl/mage/genotypes/grch38.matched"))
    ph = {r: pd.read_csv(ROOT / f"results/05-eqtl/mage/{r}/star/phenotypes.bed.gz", sep="\t", index_col=3).iloc[:, 3:]
          for r in ("chm13v2.0", "grch38", "grch38-ffx")}
    rows = []
    for sym in EQTL_EXAMPLES:
        x = g[g.symbol == sym].iloc[0]
        if sym == "PWP2":   # the FixItFelix artefact: GRCh38 lead, unpatched vs patched GRCh38
            pairs = [("grch38", x["variant_id_grch38"], gb), ("grch38-ffx", x["variant_id_grch38"], gb)]
        else:
            ee = e[e.symbol == sym].iloc[0]
            pairs = [("chm13v2.0", ee.var_a, ga), ("grch38-ffx", ee.var_b, gb)]
        for ref, vid, reader in pairs:
            geno = pd.Series(reader.read(vid), index=reader.sample_ids)
            expr = ph[ref].loc[x[f"phenotype_id_{ref}"]]
            s = expr.index.intersection(geno.index)
            rows.append(pd.DataFrame({"symbol": sym, "reference": ref, "variant": vid, "sample": s,
                                      "genotype": geno[s].values, "expression": expr[s].values}))
    ex = pd.concat(rows)
    # code genotypes as the T2T ALT dosage so both panels read the same way (sign from genotype correlation)
    for sym, d in ex.groupby("symbol"):
        refs = d.reference.unique()
        a = d[d.reference == refs[0]].set_index("sample").genotype
        b = d[d.reference == refs[1]].set_index("sample").genotype
        if np.corrcoef(a, b.reindex(a.index))[0, 1] < 0:
            ex.loc[(ex.symbol == sym) & (ex.reference == refs[1]), "genotype"] = 2 - ex.loc[(ex.symbol == sym) & (ex.reference == refs[1]), "genotype"]
    ex.to_csv(OUT / "eqtl-examples.tsv", sep="\t", index=False)
    for f in sorted(OUT.glob("*.tsv")):
        print(f.name, sum(1 for _ in open(f)) - 1, "rows")


if __name__ == "__main__":
    main()
