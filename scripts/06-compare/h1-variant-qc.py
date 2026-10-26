"""Step 6: quality check of the H1 lead variants (eGenes whose T2T lead variant lies in sequence GRCh38 lacks).

    pixi run python scripts/06-compare/h1-variant-qc.py [--cohort mage] [--method star]

The phased 1KGP T2T panel has no per-variant quality scores (only PASS sites and genotypes), so quality is judged from:
  genotypes   heterozygosity within each population (excess = reads of two paralogues called as one site),
              Hardy-Weinberg p within populations (Fisher's combination), multiallelic site in the 3202 panel
  region      censat class, segmental duplication (identity bin from the track colour), RepeatMasker class,
              variant type (SNV, indel, repeat-unit indel), inside an exon of the gene (allelic mapping bias)
  signal      best shared (matched) variant for the gene, its LD r2 with the lead, and the lead's p-value after
              conditioning on it (expression and genotypes residualised on the eQTL covariates)
Input: results/06-compare/eqtl-comparison/<cohort>/h1-t2t-only-lead-variants.tsv and the step-5 outputs.
Writes results/06-compare/eqtl-comparison/<cohort>/h1-variant-qc.tsv and prints a summary. --classify-only re-applies
the calls to that table without recomputing (the full run takes about 45 min, mostly reading genotypes).
"""
import argparse
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from tensorqtl import pgen

ROOT = Path(__file__).resolve().parents[2]
ANN = ROOT / "data/references/t2t-annotations"
PANEL = next((ROOT / "data/external/1kgp-chm13v2.0-phased-shapeit5-v1.1").glob("*.3202.bcf.gz"))


def overlaps(bed, chrom, start, end, cols):
    """Rows of a BED-like frame overlapping [start, end) on chrom."""
    b = bed[(bed[0] == chrom) & (bed[1] < end) & (bed[2] > start)]
    return b[cols]


def sd_identity(rgb):
    # T2T SD track: orange > 99%, yellow 98-99%, grey 90-98%
    return {"255,103,0": ">99%", "204,204,0": "98-99%"}.get(rgb, "90-98%")


def variant_type(ref, alt):
    if len(ref) == len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "MNV"
    short, long_ = sorted((ref, alt), key=len)
    ins = long_[len(short):] if long_.startswith(short) else long_[: len(long_) - len(short) + 1]
    # repeat-unit indel: the inserted/deleted sequence is a tandem copy of a short motif (STR / satellite)
    for k in range(1, 7):
        if len(ins) >= 2 * k and ins == (ins[:k] * (len(ins) // k + 1))[: len(ins)]:
            return f"repeat indel ({ins[:k]})n" if k > 1 else f"homopolymer indel ({ins[0]})n"
    return "indel"


def het_test(g, pops):
    """Inbreeding coefficient F = 1 - O/E over populations, and Fisher-combined within-population HWE p."""
    obs = exp = 0.0
    chi = []
    for p in pops.unique():
        x = g[pops.values == p]
        x = x[~np.isnan(x)]
        n = len(x)
        if n < 10:
            continue
        q = x.mean() / 2
        e = 2 * q * (1 - q) * n
        o = (x == 1).sum()
        obs += o; exp += e
        if 0 < q < 1:
            n0, n1, n2 = (x == 0).sum(), o, (x == 2).sum()
            ex = np.array([(1 - q) ** 2, e / n, q ** 2]) * n
            chi.append(stats.chi2.sf(((np.array([n0, n1, n2]) - ex) ** 2 / np.maximum(ex, 1e-9)).sum(), 1))
    f = 1 - obs / exp if exp else np.nan
    p = stats.combine_pvalues(chi, method="fisher")[1] if chi else np.nan
    return f, p


def residualise(y, c):
    beta, *_ = np.linalg.lstsq(c, y, rcond=None)
    return y - c @ beta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", default="mage")
    ap.add_argument("--method", default="star")
    ap.add_argument("--classify-only", action="store_true", help="re-apply the calls to an existing h1-variant-qc.tsv")
    args = ap.parse_args()
    out = ROOT / f"results/06-compare/eqtl-comparison/{args.cohort}"
    if args.classify_only:
        t = classify(pd.read_csv(out / "h1-variant-qc.tsv", sep="\t"))
        t.to_csv(out / "h1-variant-qc.tsv", sep="\t", index=False, float_format="%.4g")
        report(t)
        return
    h1 = pd.read_csv(out / "h1-t2t-only-lead-variants.tsv", sep="\t")
    d = ROOT / f"results/05-eqtl/{args.cohort}/chm13v2.0/{args.method}"

    # samples, ancestry, covariates, expression
    runs = pd.read_csv(ROOT / f"results/03-expression/{args.cohort}/matrices/runs.tsv", sep="\t")
    pop = runs.drop_duplicates("sample").set_index("sample")["ancestry"].str.split(" ").str[0]
    ph = pd.read_csv(d / "phenotypes.bed.gz", sep="\t").set_index("phenotype_id")
    cov = pd.read_csv(d / "covariates.tsv", sep="\t", index_col=0).T
    gdir = ROOT / f"results/05-eqtl/{args.cohort}/genotypes"
    rd = pgen.PgenReader(str(gdir / "t2t.all"))
    samples = [s for s in rd.sample_ids if s in cov.index and s in ph.columns]
    sidx = [rd.sample_ids.index(s) for s in samples]
    C = np.column_stack([np.ones(len(samples)), cov.loc[samples].values])
    pops = pop.reindex(samples)
    matched_ids = set(pd.read_csv(gdir / "t2t.matched.pvar", sep="\t", usecols=["ID"])["ID"])

    # annotations
    sd = pd.read_csv(ANN / "chm13v2.0_SD.bed", sep="\t", header=None, usecols=[0, 1, 2, 3, 8])
    cs = pd.read_csv(ANN / "chm13v2.0_censat_v2.1.bed", sep="\t", header=None, usecols=[0, 1, 2, 3])
    names = pd.read_csv(ROOT / "config/chrom-names-chm13v2.0.tsv", sep="\t").set_index("refseq")["ucsc"]
    gtf = pd.read_csv(ROOT / "data/star-references/chm13v2.0/annotation.gtf", sep="\t", header=None, comment="#",
                      usecols=[0, 2, 3, 4, 8])
    gtf = gtf[gtf[2] == "exon"]
    gtf["gene"] = gtf[8].str.extract(r'gene_id "([^"]+)"')[0]
    gtf[0] = gtf[0].map(names)
    gtf = gtf[gtf.gene.isin(h1.phenotype_id)]

    # RepeatMasker: only the lead-variant chromosomes, read once
    loci = h1.variant_id.str.split(":", expand=True)
    rm_rows = []
    chroms = set(loci[0])
    for chunk in pd.read_csv(ANN / "chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed", sep="\t", header=None,
                             usecols=[0, 1, 2, 3, 6, 7], chunksize=1_000_000):
        rm_rows.append(chunk[chunk[0].isin(chroms)])
    rm = pd.concat(rm_rows)

    # multiallelic in the 3202 panel
    regions = "\n".join(f"{c}\t{p}\t{p}" for c, p in zip(loci[0], loci[1]))
    rf = out / "h1-regions.tmp"
    rf.write_text(regions + "\n")
    q = subprocess.run(["bcftools", "query", "-R", str(rf), "-f", "%CHROM:%POS\t%REF\t%ALT\n", str(PANEL)],
                       capture_output=True, text=True, check=True).stdout
    rf.unlink()
    nalt = {}
    for line in q.strip().splitlines():
        k, _, alt = line.split("\t")
        nalt[k] = max(nalt.get(k, 0), len(alt.split(",")))

    nominal = {}
    rows = []
    for r in h1.itertuples():
        chrom, pos, ref, alt = r.variant_id.split(":", 3)
        pos = int(pos); end = pos + len(ref)
        rec = {"gene": r.phenotype_id, "variant_id": r.variant_id, "distance_to_tss": r.start_distance,
               "maf": min(r.af, 1 - r.af), "slope": r.slope, "pval_beta": r.pval_beta, "qval": r.qval,
               "qval_matched": r.qval_matched, "type": variant_type(ref, alt),
               "panel_alt_alleles": nalt.get(f"{chrom}:{pos}", 0)}
        c = overlaps(cs, chrom, pos, end, [3])
        rec["censat"] = ";".join(sorted(set(c[3].str.split("_").str[0]))) if len(c) else ""
        s = overlaps(sd, chrom, pos, end, [8])
        rec["sd_identity"] = min((sd_identity(x) for x in s[8]), key=[">99%", "98-99%", "90-98%"].index) \
            if len(s) else ""
        m = overlaps(rm, chrom, pos, end, [3, 6])
        rec["repeatmasker"] = ";".join(sorted(set(m[6] + ":" + m[3]))) if len(m) else ""
        e = gtf[(gtf.gene == r.phenotype_id) & (gtf[0] == chrom) & (gtf[3] <= end) & (gtf[4] >= pos)]
        rec["in_gene_exon"] = len(e) > 0

        g = np.asarray(rd.read(r.variant_id), dtype=float)[sidx]
        rec["het_F"], rec["hwe_p"] = het_test(g, pops)

        # best shared variant for this gene, from the all-variants nominal output
        if chrom not in nominal:
            nominal[chrom] = pd.read_parquet(d / f"all/nominal/{chrom}.cis_qtl_pairs.{chrom}.parquet",
                                             columns=["phenotype_id", "variant_id", "pval_nominal"],
                                             filters=[("phenotype_id", "in", list(h1.phenotype_id))])
        n = nominal[chrom]
        n = n[(n.phenotype_id == r.phenotype_id) & n.variant_id.isin(matched_ids)]
        y = residualise(ph.loc[r.phenotype_id, samples].values.astype(float), C)
        gl = residualise(np.nan_to_num(g, nan=np.nanmean(g)), C)
        if len(n):
            best = n.loc[n.pval_nominal.idxmin()]
            gs = np.asarray(rd.read(best.variant_id), dtype=float)[sidx]
            gs = residualise(np.nan_to_num(gs, nan=np.nanmean(gs)), C)
            rec["best_shared_variant"], rec["best_shared_p"] = best.variant_id, best.pval_nominal
            rec["ld_r2"] = np.corrcoef(gl, gs)[0, 1] ** 2
            if rec["ld_r2"] > 0.999:  # same genotypes: the two cannot be separated
                rows.append(rec)
                continue
            X = np.column_stack([gl, gs])
            beta, *_ = np.linalg.lstsq(X, y, rcond=None)
            res = y - X @ beta
            dof = len(y) - C.shape[1] - 2
            se = np.sqrt(res @ res / dof * np.linalg.inv(X.T @ X).diagonal())
            rec["lead_p_given_shared"] = 2 * stats.t.sf(abs(beta[0] / se[0]), dof)
            rec["shared_p_given_lead"] = 2 * stats.t.sf(abs(beta[1] / se[1]), dof)
        rows.append(rec)

    t = classify(pd.DataFrame(rows))
    t.to_csv(out / "h1-variant-qc.tsv", sep="\t", index=False, float_format="%.4g")
    report(t)


def classify(t):
    """Genotype call (is the variant reliably genotyped) and signal call (does the T2T-only variant add anything)."""
    t = t.drop(columns=["genotype_flags", "genotype_call", "signal"], errors="ignore")
    # censat classes other than ct (centromere transition) are satellite arrays
    sat = t.censat.fillna("").str.split(";").apply(lambda c: any(x not in ("", "ct") for x in c))
    # indel (or padded MNV) in a tandem repeat: by its own sequence or by RepeatMasker
    tandem = t.repeatmasker.fillna("").str.contains("Simple_repeat|Satellite|Low_complexity")
    in_rep = t.type.str.contains("repeat|homopolymer") | (~t.type.eq("SNV") & tandem)
    flags = {"het_excess": t.het_F < -0.1, "het_deficit": t.het_F > 0.1, "hwe": t.hwe_p < 1e-6,
             "multiallelic": t.panel_alt_alleles > 1, "satellite": sat, "sd>99%": t.sd_identity.eq(">99%"),
             "repeat_indel": in_rep, "exonic": t.in_gene_exon.astype(bool)}
    t["genotype_flags"] = [";".join(k for k, v in flags.items() if v.iloc[i]) for i in range(len(t))]
    # doubtful: genotypes look like merged paralogues or miscalls, or a length change inside a satellite array
    strong = flags["het_excess"] | flags["hwe"] | (sat & in_rep)
    t["genotype_call"] = np.where(strong, "doubtful", np.where(t.genotype_flags != "", "borderline", "solid"))
    lead = t.lead_p_given_shared < 0.05 / len(t)
    shared = t.shared_p_given_lead < 0.05
    t["signal"] = np.select(
        [t.best_shared_p.isna(), t.ld_r2 > 0.999, lead & ~shared, lead & shared, ~lead & shared],
        ["no shared variant", "perfect proxy shared", "T2T-only variant explains it", "both contribute",
         "explained by a shared variant"], "not separable (LD)")
    return t


def report(t):
    print(f"{len(t)} H1 lead variants ({t.variant_id.nunique()} distinct)")
    print(pd.crosstab(t.genotype_call, t.signal, margins=True).to_string())
    nm = ~(t.qval_matched <= 0.05)
    print(f"\n{nm.sum()} genes eGenes only with T2T-only variants (or not testable with matched variants):")
    print(t.loc[nm, ["gene", "variant_id", "type", "genotype_call", "genotype_flags", "het_F", "qval",
                     "qval_matched"]].to_string(index=False))

if __name__ == "__main__":
    main()
