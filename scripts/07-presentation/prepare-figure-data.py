"""Presentation: export small tables for the R figures (make-figures.R) and the workbook (build-workbook.py).

    pixi run python scripts/07-presentation/prepare-figure-data.py

Reads finished results only (MAGE discovery; Geuvadis replication). Writes results/07-presentation/data/*.tsv:
  cohorts.tsv           one row per run of both samplesheets (cohort, population, superpopulation, sex, read pairs)
  mapping-qc.tsv        MAGE runs.tsv: uniquely / multi-mapped / unmapped % per run and reference
  variant-matching.tsv  results/06-compare/variant-matching/summary.tsv
  egenes.tsv            genes tested and eGenes per reference and quantification method
  effects.tsv           same-variant comparison T2T vs patched GRCh38 (STAR), all genes, with robust flag (the 57, old
                        test) and tier (gene sets: core / main / extent, define-gene-sets.py)
  positive-control.tsv  q-values on the three references for the FixItFelix genes
  robust-genes.tsv      the robust set joined with its read destination; effect_class splits combine-methods.py's pattern:
                        opposite-sign genes are "reversed" when both slopes have p < 1e-5 (docs 13B.3), else "opposite, one weak"
  replication-h0.tsv    the 57 genes at the same variant pair in the 358 Geuvadis individuals not in MAGE
                        (replicate-changes.py --exclude-shared) with effect_class
  replication-h1.tsv    the H1 genes at the same variant in the 358 (replicate-h1.py --exclude-shared)
  gene-sets.tsv         the 255 extent genes with tier, paired-test deltas and Geuvadis replication (358 individuals)
  paired-effects.tsv    every compared gene, STAR: z_diff (independent) and z_paired, residual correlation rho, tier
  calibration.tsv       bootstrap and genomic-inflation checks of the paired test (paired-test.py)
  size-floor.tsv        gene count and Geuvadis replication for each definition of "changed" (docs 13B.4)
  replication-paired.tsv  the 255 at MAGE's variant pairs in the 358, paired test (paired-test.py --at mage)
  copy-number.tsv       extent genes with ctyper copy number: delta without / with copy number as covariate
  coloc.tsv             gene x trait colocalisation on both references (coloc.py), rows with PP4 >= 0.5 on either
  gwas-lookup.tsv       GWAS Catalog proxy hits per gene: main set and matched background (gwas-lookup.py)
  paralog-deltas.tsv    per-individual CPM change (T2T - patched GRCh38) for example gene pairs
  eqtl-examples.tsv     normalised expression and lead-variant genotype on both references for example genes
  expression-compare.tsv  mean CPM per gene on each reference (MAGE, 731, STAR and salmon) for the three pairs, with
                        region class, "changed" (>= 2-fold under STAR and salmon, same direction) and the eQTL result
                        (compared with the paired test; in the main changed-eQTL set)
  variant-matching-mage.tsv  variant pairing as seen by the MAGE eQTL scan (chr1-22, MAF >= 1% in MAGE), with SNP/indel
  expression-eqtl-examples.tsv  per-person CPM and genotype on T2T and patched GRCh38 for EXPR_EXAMPLES
  paired-example.tsv    PAIRED_EXAMPLE, STAR, T2T vs patched GRCh38: per-person expression after the covariates on each
                        reference (as paired-test.py) and the T2T genotype, to show why the paired test works
  top10-changed.tsv     the 10 main-set genes with the largest change (T2T vs patched GRCh38, STAR): variant pair with
                        positions, dbSNP rsID of the GRCh38 variant (Ensembl REST, GRCh38; cached in
                        data/metadata/rsid-cache.tsv), GRCh38 REF/ALT, the effect allele (the one the slopes refer
                        to) and its frequency in MAGE, 1000 Genomes (ALL) and gnomAD genomes (cached in
                        data/metadata/af-cache.tsv), slopes, change, paired z and log10 p, Geuvadis replication
  changed-mechanism.tsv  the main-set genes (T2T vs patched GRCh38) by likely reason for the change, first that applies:
                        sd98 (>= 98% segmental duplication) > one_ref_partner (read partner annotated on one reference
                        only, r <= -0.5) > sd90 (90-98% duplication covering >= 50% of the gene) > paralog_partner (read
                        partner r <= -0.7) > traced_no_partner (read destination run, no clear partner) > not_traced
                        (read destination not yet run: the genes added by the paired test)
  coloc-table.tsv       main-set genes that colocalise on one reference only (coloc.tsv): the top shared variant on the
                        reference where they colocalise (GRCh38 ID, as coloc.py tests both references on GRCh38 IDs),
                        its dbSNP rsID, and GWAS P at that variant
  coverage-example.tsv  per-base RNA coverage (uniquely mapped reads, per million) along the exons of COVERAGE_PAIR on
                        T2T and patched GRCh38, mean per genotype group of the first gene's lead variant, from the ASE
                        target BAMs (results/06-compare/ase/rna/<run>/<ref>/targets.bam; samtools depth -a -Q 255 -r)
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
EXPR_EXAMPLES = ["PKD1P4-NPIPA8", "CCZ1", "ASAH2"]   # expression changes, eQTL does not / the reverse; ASAH2: Fig 22
COVERAGE_PAIR = ("ASAH2", "ASAH2B")   # paralog pair whose reads move between references
PAIRED_EXAMPLE = "IFITM2"   # main-set gene, rho 0.99: missed by the independent test (z -3.1), paired z -25


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

    pd.read_csv(ROOT / "results/04-expression/mage/matrices/runs.tsv", sep="\t").to_csv(OUT / "mapping-qc.tsv", sep="\t", index=False)
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
    gs = pd.read_csv(CMP / "gene-sets/changed-chm13v2.0-vs-grch38-ffx.tsv", sep="\t", dtype={"gene_ncbi": str})
    e = e.merge(gs[["gene_ncbi", "tier"]], on="gene_ncbi", how="left")
    e["tier"] = e.tier.fillna("unchanged")
    e.to_csv(OUT / "effects.tsv", sep="\t", index=False)

    g = pd.read_csv(CMP / "star-matched/genes.tsv", sep="\t", dtype={"gene_ncbi": str}, low_memory=False)
    pc = g[g.symbol.isin(["PWP2", "GATD3", "CBS", "U2AF1", "ICOSLG"])]
    pc = pc.melt(id_vars=["symbol"], value_vars=[f"qval_{r}" for r in REFS], var_name="reference", value_name="qval")
    pc["reference"] = pc.reference.str.replace("qval_", "")
    pc.to_csv(OUT / "positive-control.tsv", sep="\t", index=False)

    rd = pd.read_csv(CMP / "robust-matched/read-destination-chm13v2.0-vs-grch38-ffx.tsv", sep="\t", dtype={"gene_ncbi": str})
    rg = robust[~robust.artefact_flag].merge(rd.drop(columns=["symbol", "region_class", "stronger_on"]), on="gene_ncbi", how="left")
    strong = (rg.log10p_a > 5) & (rg.log10p_b > 5)
    rg["effect_class"] = np.select(
        [rg.pattern.eq("opposite sign") & strong, rg.pattern.eq("opposite sign"), rg.pattern.eq("stronger on chm13v2.0")],
        ["Reversed", "Opposite sign, one weak", "Stronger on T2T"], "Stronger on patched GRCh38")
    rg.to_csv(OUT / "robust-genes.tsv", sep="\t", index=False)
    rep = CMP.parent / "replication-geuvadis-independent"
    h0 = pd.read_csv(rep / "genes-chm13v2.0-vs-grch38-ffx.tsv", sep="\t", dtype={"gene_ncbi": str})
    h0 = h0.merge(rg[["gene_ncbi", "effect_class", "pattern"]], on="gene_ncbi", how="left")
    h0.to_csv(OUT / "replication-h0.tsv", sep="\t", index=False)
    pd.read_csv(rep / "h1-genes.tsv", sep="\t").to_csv(OUT / "replication-h1.tsv", sep="\t", index=False)

    paired_tables(gs)

    # per-individual CPM deltas for example pairs (same construction as read-destination.py)
    cpm = {}
    for r in ("chm13v2.0", "grch38-ffx"):
        m = pd.read_csv(ROOT / f"results/04-expression/mage/matrices/{r}/star-counts.tsv.gz", sep="\t", index_col=0)
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
    expression_tables()
    variant_matching_mage()
    paired_example()
    top10_table()
    changed_mechanism()
    coloc_table()
    coverage_example()
    for f in sorted(OUT.glob("*.tsv")):
        print(f.name, sum(1 for _ in open(f)) - 1, "rows")


def paired_tables(gs):
    """Tables for the paired test, gene sets, copy number and GWAS (docs 13B.3, 13B.4)."""
    gs.to_csv(OUT / "gene-sets.tsv", sep="\t", index=False)
    pe = pd.read_csv(CMP / "paired-matched/effects-chm13v2.0-vs-grch38-ffx-star.tsv", sep="\t", dtype={"gene_ncbi": str})
    pe = pe.merge(gs[["gene_ncbi", "tier"]], on="gene_ncbi", how="left")
    pe["tier"] = pe.tier.fillna("unchanged")
    pe[["gene_ncbi", "symbol", "region_class", "artefact_flag", "slope_a", "slope_b", "rho", "z_diff", "z_paired",
        "tier"]].to_csv(OUT / "paired-effects.tsv", sep="\t", index=False)
    pd.read_csv(CMP / "paired-matched/calibration.tsv", sep="\t").to_csv(OUT / "calibration.tsv", sep="\t", index=False)

    # replication by the size floor on |delta| (decided 2026-10-29: 0.1 under both methods)
    d, ds = gs.delta.abs(), gs.delta_salmon.abs()
    rows = []
    for name, keep in [("Paired, no floor (extent)", d >= 0), ("Paired, |delta| >= 0.05 (STAR)", d >= 0.05),
                       ("Paired, |delta| >= 0.1 (STAR)", d >= 0.1), ("Paired, |delta| >= 0.1 both (main)", (d >= 0.1) & (ds >= 0.1)),
                       ("Paired, |delta| >= 0.2 both", (d >= 0.2) & (ds >= 0.2)), ("Old test (core)", gs.tier.eq("core"))]:
        x = gs[keep]
        t = x.geuvadis_same_sign.notna()
        rows.append({"definition": name, "genes": len(x), "of_57": int(x.tier.eq("core").sum()),
                     "sd": int(x.region_class.eq("sd").sum()), "tested": int(t.sum()),
                     "same_sign": int(x.geuvadis_same_sign.eq(True).sum()), "nominal": int(x.geuvadis_nominal.eq(True).sum()),
                     "bonferroni": int(x.geuvadis_bonferroni.eq(True).sum())})
    pd.DataFrame(rows).to_csv(OUT / "size-floor.tsv", sep="\t", index=False)

    rp = pd.read_csv(CMP.parent / "geuvadis/paired-matched-at-mage/replication-chm13v2.0-vs-grch38-ffx.tsv", sep="\t",
                     dtype={"gene_ncbi": str})
    rp = rp.merge(gs[["gene_ncbi", "tier", "z_paired"]].rename(columns={"z_paired": "mage_z_paired"}), on="gene_ncbi")
    rp.to_csv(OUT / "replication-paired.tsv", sep="\t", index=False)

    cn = pd.read_csv(CMP / "paired-matched-cn/effects-chm13v2.0-vs-grch38-ffx-star.tsv", sep="\t", dtype={"gene_ncbi": str})
    cr = pd.read_csv(CMP / "paired-matched-cn/robust-chm13v2.0-vs-grch38-ffx.tsv", sep="\t", dtype={"gene_ncbi": str})
    c = gs[gs.gene_ncbi.isin(cn.gene_ncbi)][["gene_ncbi", "symbol", "tier", "region_class", "delta", "z_paired"]]
    c = c.merge(cn[["gene_ncbi", "delta", "z_paired", "cn_r2"]], on="gene_ncbi", suffixes=("", "_cn"))
    c = c.merge(cr[["gene_ncbi", "robust_paired"]], on="gene_ncbi", how="left")
    c["still_changed"] = c.robust_paired.eq(True) & (np.sign(c.delta) == np.sign(c.delta_cn))
    c.drop(columns="robust_paired").to_csv(OUT / "copy-number.tsv", sep="\t", index=False)

    gw = ROOT / "results/06-compare/gwas"
    co = pd.read_csv(gw / "coloc-mage-chm13v2.0-vs-grch38-ffx.tsv", sep="\t", dtype={"gene_ncbi": str})
    co[co[["PP4_chm13v2.0", "PP4_grch38-ffx"]].max(axis=1) >= 0.5].to_csv(OUT / "coloc.tsv", sep="\t", index=False)
    pd.concat([pd.read_csv(gw / "genes-mage-chm13v2.0-vs-grch38-ffx.tsv", sep="\t"),
               pd.read_csv(gw / "genes-mage-chm13v2.0-vs-grch38-ffx-background.tsv", sep="\t")]
              ).to_csv(OUT / "gwas-lookup.tsv", sep="\t", index=False)


def expression_tables():
    """Expression level on each reference, before any eQTL (docs 12, step 7).

    Per reference and method: CPM = counts / the individual's total counts x 1e6, then the mean over the 731 individuals.
    Genes: one symbol per GeneID, annotated on both assemblies, no MHC/PAR/chrM/rRNA flag. Expressed in a pair: mean
    CPM >= 1 on either reference (STAR). Changed: |log2 ratio| >= 1 (2-fold, pseudocount 0.5 CPM) under STAR and salmon,
    same direction."""
    gc = pd.read_csv(ROOT / "results/06-compare/gene-classes.tsv", sep="\t", dtype={"gene_ncbi": str})
    gc = gc[~gc.symbol.duplicated(keep=False) & gc.annotated_on.eq("both")
            & ~gc[["flag_mhc", "flag_par", "flag_chrM", "flag_rrna"]].any(axis=1)]
    mean = {}
    for r in REFS:
        for m in ("star", "salmon"):
            x = pd.read_csv(ROOT / f"results/04-expression/mage/matrices/{r}/{m}-counts.tsv.gz", sep="\t", index_col=0)
            mean[r, m] = (x / x.sum() * 1e6).mean(axis=1)
    rows = []
    for a, b in [("chm13v2.0", "grch38"), ("grch38-ffx", "grch38"), ("chm13v2.0", "grch38-ffx")]:
        pair = f"{a}-vs-{b}"
        d = gc[["gene_ncbi", "symbol", "biotype", "region_class"]].copy()
        for m in ("star", "salmon"):
            d[f"cpm_a_{m}"] = d.symbol.map(mean[a, m]); d[f"cpm_b_{m}"] = d.symbol.map(mean[b, m])
            d[f"log2fc_{m}"] = np.log2((d[f"cpm_a_{m}"] + 0.5) / (d[f"cpm_b_{m}"] + 0.5))
        d = d.dropna()
        d = d[d[["cpm_a_star", "cpm_b_star"]].max(axis=1) >= 1].copy()
        d["changed"] = (d.log2fc_star.abs() >= 1) & (d.log2fc_salmon.abs() >= 1) & (np.sign(d.log2fc_star) == np.sign(d.log2fc_salmon))
        cmp = pd.read_csv(CMP / f"paired-matched/effects-{pair}-star.tsv", sep="\t", usecols=["gene_ncbi"], dtype=str)
        gs = pd.read_csv(CMP / f"gene-sets/changed-{pair}.tsv", sep="\t", dtype={"gene_ncbi": str})
        d["eqtl_compared"] = d.gene_ncbi.isin(cmp.gene_ncbi)
        d["eqtl_changed"] = d.gene_ncbi.isin(gs.loc[gs.tier.isin(["core", "main"]), "gene_ncbi"])
        d.insert(0, "pair", pair)
        rows.append(d)
    pd.concat(rows).to_csv(OUT / "expression-compare.tsv", sep="\t", index=False, float_format="%.4g")

    # two genes for "expression change is not eQTL change" (T2T vs patched GRCh38, STAR): per-person CPM by genotype of
    # the compared variant pair (T2T ALT dosage on both panels, sign from genotype correlation)
    from tensorqtl import pgen
    e = pd.read_csv(CMP / "paired-matched/effects-chm13v2.0-vs-grch38-ffx-star.tsv", sep="\t")
    rd = {"chm13v2.0": pgen.PgenReader(str(ROOT / "results/05-eqtl/mage/genotypes/t2t.matched")),
          "grch38-ffx": pgen.PgenReader(str(ROOT / "results/05-eqtl/mage/genotypes/grch38.matched"))}
    cpm = {}
    for r in rd:
        x = pd.read_csv(ROOT / f"results/04-expression/mage/matrices/{r}/star-counts.tsv.gz", sep="\t", index_col=0)
        cpm[r] = x.loc[EXPR_EXAMPLES] / x.sum() * 1e6
    rows = []
    for sym in EXPR_EXAMPLES:
        x = e[e.symbol == sym].iloc[0]
        g = {r: pd.Series(rd[r].read(v), index=rd[r].sample_ids) for r, v in (("chm13v2.0", x.var_a), ("grch38-ffx", x.var_b))}
        ids = cpm["chm13v2.0"].columns.intersection(g["chm13v2.0"].index)
        if np.corrcoef(g["chm13v2.0"][ids], g["grch38-ffx"][ids])[0, 1] < 0:
            g["grch38-ffx"] = 2 - g["grch38-ffx"]
        for r in rd:
            rows.append(pd.DataFrame({"symbol": sym, "reference": r, "variant": x.var_a, "sample": ids,
                                      "genotype": g[r][ids].values, "cpm": cpm[r].loc[sym, ids].values,
                                      "slope": x.slope_a if r == "chm13v2.0" else x.slope_b}))
    pd.concat(rows).to_csv(OUT / "expression-eqtl-examples.tsv", sep="\t", index=False, float_format="%.4g")


def variant_matching_mage():
    """Variant pairing as seen by the MAGE eQTL scan (chr1-22, in-cohort MAF >= 0.01; docs 13B.1).

    Pairing (match-variants.sbatch) ran once on all 1,089 people with RNA-seq; MAGE's genotype files then keep the
    variants with MAF >= 0.01 within its 731 people. Each T2T variant in t2t.all gets the pairing category of the
    1,089-person table; variants below 1% in the 1,089 were never paired. GRCh38 variants without a matched T2T partner
    = grch38.all minus grch38.matched. SNP = single-base REF and ALT."""
    g = ROOT / "results/05-eqtl/mage/genotypes"
    rd = lambda f: pd.read_csv(g / f"{f}.pvar", sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4],
                               names=["chrom", "id", "ref", "alt"])
    auto = [f"chr{i}" for i in range(1, 23)]
    ta, tm, ga, gm = (x[x.chrom.isin(auto)] for x in (rd("t2t.all"), rd("t2t.matched"), rd("grch38.all"), rd("grch38.matched")))
    cat = pd.concat(pd.read_csv(ROOT / f"results/06-compare/variant-matching/{c}.tsv.gz", sep="\t",
                                usecols=["side", "t2t_id", "category"], low_memory=False).query("side == 't2t'")
                    for c in auto).set_index("t2t_id").category
    ta = ta.assign(category=ta.id.map(cat).fillna("not_paired_rare_in_1089"), side="t2t")
    ta.loc[ta.id.isin(tm.id), "category"] = "matched"
    gx = ga[~ga.id.isin(gm.id)].assign(category="grch38_no_matched_partner", side="grch38")
    d = pd.concat([ta, gx])
    d["type"] = ((d.ref.str.len() == 1) & (d.alt.str.len() == 1)).map({True: "snp", False: "indel"})
    out = d.groupby(["side", "category", "type"]).size().rename("n").reset_index()
    tot = pd.DataFrame([{"side": n, "category": "all", "type": t, "n": int(((x.ref.str.len() == 1) & (x.alt.str.len() == 1)).eq(t == "snp").sum())}
                        for n, x in (("t2t", ta), ("grch38", ga)) for t in ("snp", "indel")])
    pd.concat([out, tot]).to_csv(OUT / "variant-matching-mage.tsv", sep="\t", index=False)


def paired_example(a="chm13v2.0", b="grch38-ffx", m="star"):
    """Per-person data behind one paired test (PAIRED_EXAMPLE), computed as paired-test.py does."""
    from tensorqtl import pgen
    e = pd.read_csv(CMP / f"paired-matched/effects-{a}-vs-{b}-{m}.tsv", sep="\t", dtype={"gene_ncbi": str})
    r = e[e.symbol == PAIRED_EXAMPLE].iloc[0]
    gid = {x: pd.read_csv(ROOT / f"data/star-references/{x}/genes.tsv", sep="\t", dtype={"gene_ncbi": str})
           .drop_duplicates("gene_ncbi").set_index("gene_ncbi").gene_id for x in (a, b)}
    ph = {x: pd.read_csv(ROOT / f"results/05-eqtl/mage/{x}/{m}/phenotypes.bed.gz", sep="\t")
          .set_index("phenotype_id").iloc[:, 3:] for x in (a, b)}
    cov = {x: pd.read_csv(ROOT / f"results/05-eqtl/mage/{x}/{m}/covariates.tsv", sep="\t", index_col=0).T for x in (a, b)}
    rd = pgen.PgenReader(str(ROOT / "results/05-eqtl/mage/genotypes/t2t.matched"))
    s = [x for x in ph[a].columns if x in ph[b].columns and x in cov[a].index and x in cov[b].index and x in rd.sample_ids]
    g = pd.DataFrame(rd.read_list([r.var_a], impute_mean=True), columns=rd.sample_ids).iloc[0][s].values
    y = {}
    for x in (a, b):
        Q, _ = np.linalg.qr(np.column_stack([np.ones(len(s)), cov[x].loc[s].values]))
        v = ph[x].loc[gid[x][r.gene_ncbi], s].values.astype(float)
        y[x] = v - Q @ (Q.T @ v)
    pd.DataFrame({"sample": s, "genotype": g, "expr_t2t": y[a], "expr_patched": y[b],
                  "symbol": PAIRED_EXAMPLE, "slope_t2t": r.slope_a, "slope_patched": r.slope_b, "rho": r.rho,
                  "se_indep": r.se_diff_indep, "se_paired": r.se_diff_paired, "z_indep": r.z_diff,
                  "z_paired": r.z_paired}).to_csv(OUT / "paired-example.tsv", sep="\t", index=False, float_format="%.5g")


def rsid_lookup(variants):
    """dbSNP rsIDs for GRCh38 variant IDs chr:pos:ref:alt (Ensembl REST overlap, allele-matched); cached on disk."""
    import json, time, urllib.request
    cache_f = ROOT / "data/metadata/rsid-cache.tsv"
    cache = pd.read_csv(cache_f, sep="\t", index_col=0).rsid.to_dict() if cache_f.exists() else {}
    for v in variants:
        if v in cache:
            continue
        c, pos, ref, alt = v.split(":")
        pos = int(pos)
        end = pos + max(len(ref), len(alt)) + 20   # indels may be shifted right in repeats by Ensembl
        url = f"https://rest.ensembl.org/overlap/region/human/{c[3:]}:{pos}-{end}?feature=variation"
        hit = ""
        for _ in range(3):
            try:
                req = urllib.request.Request(url, headers={"Content-Type": "application/json", "User-Agent": "t2t-eqtl"})
                rows = json.load(urllib.request.urlopen(req, timeout=30))
                break
            except Exception:
                rows = []; time.sleep(2)
        snv = len(ref) == len(alt) == 1
        for r in rows:
            al = set(r.get("alleles", []))
            if snv and r["start"] == pos and ref in al and alt in al:
                hit = r["id"]; break
            if not snv and "-" in al:   # indel: same inserted/deleted sequence, allowing a shift within a repeat
                seq = alt[len(ref):] if len(alt) > len(ref) else ref[len(alt):]
                if any(len(a) == len(seq) and a in seq + seq for a in al - {"-"}):
                    hit = r["id"]; break
        cache[v] = hit
        time.sleep(0.1)
    pd.Series(cache, name="rsid").rename_axis("variant").to_csv(cache_f, sep="\t")
    return {v: cache[v] for v in variants}


def global_af(rsids):
    """1000 Genomes phase 3 (ALL, 2,504 people) and gnomAD genomes allele frequencies per rsID (Ensembl REST);
    cached in data/metadata/af-cache.tsv as rsid, source, allele, frequency."""
    import json, time, urllib.request
    cache_f = ROOT / "data/metadata/af-cache.tsv"
    cache = pd.read_csv(cache_f, sep="\t", keep_default_na=False) if cache_f.exists() else \
        pd.DataFrame(columns=["rsid", "source", "allele", "frequency"])
    new = []
    for rs in rsids:
        if not rs or rs in set(cache.rsid):
            continue
        for _ in range(3):
            try:
                req = urllib.request.Request(f"https://rest.ensembl.org/variation/human/{rs}?pops=1",
                                             headers={"Content-Type": "application/json", "User-Agent": "t2t-eqtl"})
                d = json.load(urllib.request.urlopen(req, timeout=30)); break
            except Exception:
                d = {}; time.sleep(2)
        for q in d.get("populations", []):
            if q["population"] in ("1000GENOMES:phase_3:ALL", "gnomADg:ALL"):
                new.append({"rsid": rs, "source": q["population"], "allele": q["allele"], "frequency": q["frequency"]})
        time.sleep(0.1)
    cache = pd.concat([cache, pd.DataFrame(new)], ignore_index=True)
    cache.to_csv(cache_f, sep="\t", index=False)
    return cache


def top10_table(n=10):
    """The n main-set genes with the largest change (T2T vs patched GRCh38, STAR), with their variant and statistics."""
    from scipy.stats import norm
    g = pd.read_csv(CMP / "gene-sets/changed-chm13v2.0-vs-grch38-ffx.tsv", sep="\t")
    g = g[g.tier.isin(["core", "main"])].assign(ad=lambda d: d.delta.abs()).sort_values("ad", ascending=False).head(n)
    rs = rsid_lookup(g.var_b.tolist())
    g["rsid"] = g.var_b.map(rs)
    g["pos_t2t"] = g.var_a.str.split(":").str[:2].str.join(":")
    g["pos_grch38"] = g.var_b.str.split(":").str[:2].str.join(":")
    g["log10p_paired"] = (norm.logsf(g.z_paired.abs()) + np.log(2)) / np.log(10)
    # alleles on GRCh38 (dbSNP orientation) and the effect allele: the allele the slopes refer to (T2T ALT), written in
    # GRCh38 terms via the genotype correlation (geno_r > 0: GRCh38 ALT; < 0: GRCh38 REF)
    al = g.var_b.str.split(":")
    g["ref_grch38"], g["alt_grch38"] = al.str[2], al.str[3]
    g["effect_allele"] = np.where(g.geno_r > 0, g.alt_grch38, g.ref_grch38)
    from tensorqtl import pgen
    rd = pgen.PgenReader(str(ROOT / "results/05-eqtl/mage/genotypes/grch38.matched"))
    alt_f = rd.read_list(g.var_b.tolist(), impute_mean=True).mean(axis=1).values / 2
    g["eaf_mage"] = np.where(g.geno_r > 0, alt_f, 1 - alt_f)
    af = global_af(g.rsid.tolist())
    def ens_allele(r):   # Ensembl writes indels without the VCF anchor base: '-' for the short allele
        ref, alt, e = r.ref_grch38, r.alt_grch38, r.effect_allele
        if len(ref) == len(alt) == 1:
            return e
        long_, short = (alt, ref) if len(alt) > len(ref) else (ref, alt)
        return "-" if e == short else long_[len(short):]
    def lookup(r, src):
        x = af[(af.rsid == r.rsid) & (af.source == src)]
        a = ens_allele(r)
        hit = x[x.allele == a]
        if hit.empty and a != "-":   # indel written rotated in a repeat
            hit = x[(x.allele != "-") & x.allele.map(lambda y: len(y) == len(a) and y in a + a)]
        return float(hit.frequency.iloc[0]) if len(hit) else np.nan
    g["eaf_1kg"] = [lookup(r, "1000GENOMES:phase_3:ALL") for r in g.itertuples()]
    g["eaf_gnomad"] = [lookup(r, "gnomADg:ALL") for r in g.itertuples()]
    g[["symbol", "biotype", "region_class", "var_a", "var_b", "pos_t2t", "pos_grch38", "rsid", "ref_grch38", "alt_grch38",
       "effect_allele", "eaf_mage", "eaf_1kg", "eaf_gnomad", "geno_r", "slope_a", "slope_b",
       "delta", "z_paired", "log10p_paired", "delta_salmon", "z_paired_salmon", "geuvadis_nominal"]].to_csv(
        OUT / "top10-changed.tsv", sep="\t", index=False, float_format="%.4g")


def changed_mechanism():
    """Main-set genes by likely reason for the eQTL change (duplication map for all; read destination for the core)."""
    g = pd.read_csv(CMP / "gene-sets/changed-chm13v2.0-vs-grch38-ffx.tsv", sep="\t")
    g = g[g.tier.isin(["core", "main"])]
    c = pd.read_csv(ROOT / "results/06-compare/gene-classes.tsv", sep="\t")[["gene_ncbi", "t2t_sd_frac"]]
    rd = pd.read_csv(CMP / "robust-matched/read-destination-chm13v2.0-vs-grch38-ffx.tsv", sep="\t")[
        ["gene_ncbi", "partner", "partner_annotated_on", "partner_r"]]
    m = g.merge(c, on="gene_ncbi", how="left").merge(rd, on="gene_ncbi", how="left")
    def cls(r):
        if r.region_class == "sd":
            return "sd98"
        if r.partner_annotated_on in ("t2t_only", "grch38_only") and r.partner_r <= -0.5:
            return "one_ref_partner"
        if r.t2t_sd_frac >= 0.5:
            return "sd90"
        if r.partner_r <= -0.7:
            return "paralog_partner"
        return "traced_no_partner" if pd.notna(r.partner) else "not_traced"
    m["mechanism"] = m.apply(cls, axis=1)
    m[["gene_ncbi", "symbol", "tier", "region_class", "t2t_sd_frac", "partner", "partner_annotated_on", "partner_r",
       "mechanism"]].to_csv(OUT / "changed-mechanism.tsv", sep="\t", index=False, float_format="%.4g")


def coloc_table():
    """Main-set genes colocalising on one reference only, with the top shared variant and its GWAS P."""
    from scipy.stats import norm
    c = pd.read_csv(OUT / "coloc.tsv", sep="\t")
    c = c[c.tier.isin(["core", "main"]) & c.call.str.endswith("only", na=False)].copy()
    c["on"] = np.where(c.call == "chm13v2.0 only", "chm13v2.0", "grch38-ffx")
    c["top"] = [r[f"top_{r.on}"] for _, r in c.iterrows()]
    c["gwas_z"] = [r[f"gwas_z_{r.on}"] for _, r in c.iterrows()]
    c["gwas_log10p"] = (norm.logsf(c.gwas_z.abs()) + np.log(2)) / np.log(10)
    c["rsid"] = c.top.map(rsid_lookup(c.top.tolist()))
    c[["symbol", "trait", "gcst", "on", "top", "rsid", "gwas_z", "gwas_log10p", "PP4_chm13v2.0", "PP4_grch38-ffx"]].to_csv(
        OUT / "coloc-table.tsv", sep="\t", index=False, float_format="%.4g")


def coverage_example(pair=COVERAGE_PAIR, refs=("chm13v2.0", "grch38-ffx"), binsize=10, out_name="coverage-example.tsv"):
    """Mean per-base coverage per genotype group along the exons of a paralog pair, from the ASE target BAMs."""
    import io, shutil, subprocess
    from tensorqtl import pgen
    sam = shutil.which("samtools") or str(ROOT / ".pixi/envs/bio/bin/samtools")  # pixi environment `bio`
    ase = ROOT / "results/06-compare/ase"
    runs = pd.read_csv(ase / "targets/runs.tsv", sep="\t")
    runs = runs[[(ase / "rna" / r / "complete.txt").exists() for r in runs.run]]
    g = pd.read_csv(CMP / "gene-sets/changed-chm13v2.0-vs-grch38-ffx.tsv", sep="\t").set_index("symbol").loc[pair[0]]
    rd = pgen.PgenReader(str(ROOT / "results/05-eqtl/mage/genotypes/t2t.matched"))
    geno = pd.DataFrame(rd.read_list([g.var_a], impute_mean=True), columns=rd.sample_ids).iloc[0].round()
    rows = []
    for ref in refs:
        gtf = ROOT / f"data/star-references/{ref}/annotation.gtf"
        ex = []
        for line in open(gtf):
            f = line.split("\t")
            if len(f) > 8 and f[2] == "exon":
                for sym in pair:
                    if f'gene_id "{sym}";' in f[8]:
                        ex.append((sym, f[0], int(f[3]), int(f[4]), f[6]))
        ex = pd.DataFrame(ex, columns=["gene", "chrom", "start", "end", "strand"])
        # union of exons per gene, then concatenated exonic positions in transcript direction (5' to 3')
        segs = {}
        for sym, d in ex.groupby("gene"):
            iv = []
            for a, b in sorted(zip(d.start, d.end)):
                if iv and a <= iv[-1][1] + 1:
                    iv[-1][1] = max(iv[-1][1], b)
                else:
                    iv.append([a, b])
            segs[sym] = (d.chrom.iloc[0], d.strand.iloc[0], iv)
        spans = [f"{c}:{iv[0][0]}-{iv[-1][1]}" for c, st, iv in segs.values()]   # indexed queries (fast)
        pos_index = {}
        for sym, (c, st, iv) in segs.items():
            bases = [x for a, b in iv for x in range(a, b + 1)]
            if st == "-":
                bases = bases[::-1]
            for k, x in enumerate(bases):
                pos_index[(c, x)] = (sym, k)
        per_sample = {}
        for r in runs.itertuples():
            d = ase / "rna" / r.run / ref
            uniq = next(int(l.split("|")[1]) for l in open(d / "Log.final.out") if "Uniquely mapped reads number" in l)
            out = "".join(subprocess.run([sam, "depth", "-a", "-Q", "255", "-r", sp_, str(d / "targets.bam")],
                                         check=True, capture_output=True, text=True).stdout for sp_ in spans)
            dep = pd.read_csv(io.StringIO(out), sep="\t", header=None, names=["chrom", "pos", "depth"])
            key = [pos_index.get((c, x)) for c, x in zip(dep.chrom, dep.pos)]
            dep["gene"] = [k[0] if k else None for k in key]
            dep["i"] = [k[1] if k else -1 for k in key]
            dep = dep[dep.gene.notna()]
            dep["cpm"] = dep.depth / uniq * 1e6
            per_sample.setdefault(r.sample, []).append(dep.set_index(["gene", "i"]).cpm)
        prof = pd.DataFrame({smp: pd.concat(v, axis=1).mean(axis=1) for smp, v in per_sample.items()})
        prof = prof[[c for c in prof.columns if c in geno.index]]
        gg = geno[prof.columns]
        for k in sorted(gg.unique()):
            m = prof.loc[:, gg == k].mean(axis=1).rename("cpm").reset_index()
            m["bin"] = m.i // binsize
            m = m.groupby(["gene", "bin"], as_index=False).cpm.mean()
            m["ref"], m["genotype"], m["n"] = ref, int(k), int((gg == k).sum())
            rows.append(m)
    out = pd.concat(rows)
    out["pos_kb"] = out["bin"] * binsize / 1000
    out["variant_t2t"] = g.var_a
    out.to_csv(OUT / out_name, sep="\t", index=False, float_format="%.5g")


if __name__ == "__main__":
    main()
