"""Step 6 (H6): SuSiE-coloc at the loci where coloc.abf differs between references, allowing several causal variants.

    pixi run python scripts/06-compare/coloc-susie.py [--cohort mage] [--pair chm13v2.0 grch38-ffx] [--min-pp4 0.5]

coloc.abf (coloc.py) assumes one causal variant per trait. A colocalisation lost on one reference could instead be a
second eQTL signal; SuSiE fine-mapping (Wang et al. 2020) with coloc.susie (Wallace 2021) tests every pair of credible
sets. Loci: gene x trait rows of coloc-<cohort>-<A>-vs-<B>.tsv with PP4 >= --min-pp4 on either reference.
Per locus, +- 500 kb around the coloc.abf top variants, biallelic SNPs present in the eQTL results of both references,
the GWAS, the cohort's GRCh38 genotypes and the 1KGP EUR panel, with MAF >= 0.01 in the cohort and in EUR:
  eQTL (each reference): beta = slope, varbeta = se^2, N = cohort size, sdY = 1; LD = genotype correlation in the cohort
  GWAS: beta, se from the harmonised file; type cc, N from its metadata; LD from the 525 unrelated 1KGP EUR
  (NYGC GRCh38 panel, 20130606_g1k_3202_samples_ped_population.txt)
All on the GRCh38 alt allele of the variant ID (T2T slopes via the matching table; SNPs only).
R: coloc::runsusie per dataset, coloc.susie per reference (coloc 5; default priors), coloc.abf on the same variants.
Writes results/06-compare/gwas/coloc-susie-<cohort>-<A>-vs-<B>.tsv (max PP4 over credible-set pairs per reference).
"""
import argparse
import glob
import importlib.util
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pgenlib
import yaml

ROOT = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location("coloc", ROOT / "scripts/06-compare/coloc.py")
coloc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(coloc)
FLANK = 500_000
PANEL = ROOT / "data/external/1kgp-grch38-nygc-phased-20220422"

R_SCRIPT = r"""
suppressMessages({library(coloc); library(data.table)})
args <- commandArgs(TRUE); dir <- args[1]
loci <- fread(file.path(dir, "loci.tsv"))
out <- list()
for (i in seq_len(nrow(loci))) {
  L <- loci[i]; d <- file.path(dir, L$locus)
  v <- fread(file.path(d, "variants.tsv"))
  ldq <- as.matrix(fread(file.path(d, "ld-cohort.tsv"), header = FALSE)); ldg <- as.matrix(fread(file.path(d, "ld-eur.tsv"), header = FALSE))
  dimnames(ldq) <- list(v$snp, v$snp); dimnames(ldg) <- list(v$snp, v$snp)
  g <- list(beta = v$g_beta, varbeta = v$g_se^2, snp = v$snp, position = v$pos, type = "cc", N = L$gwas_n,
            MAF = v$maf_eur, LD = ldg)
  sg <- tryCatch(runsusie(g, repeat_until_convergence = FALSE, maxit = 1000), error = function(e) NULL)
  for (ref in strsplit(L$refs, ",")[[1]]) {
    q <- list(beta = v[[paste0("slope_", ref)]], varbeta = v[[paste0("se_", ref)]]^2, snp = v$snp, position = v$pos,
              type = "quant", N = L$eqtl_n, sdY = 1, MAF = v$maf_cohort, LD = ldq)
    sq <- tryCatch(runsusie(q, repeat_until_convergence = FALSE, maxit = 1000), error = function(e) NULL)
    abf <- coloc.abf(q[c("beta","varbeta","snp","type","N","sdY","MAF")], g[c("beta","varbeta","snp","type","N","MAF")])$summary
    rec <- data.table(locus = L$locus, ref = ref, nsnps = nrow(v), abf_PP4 = abf[["PP.H4.abf"]],
                      abf_PP3 = abf[["PP.H3.abf"]], cs_eqtl = if (is.null(sq)) NA else length(sq$sets$cs),
                      cs_gwas = if (is.null(sg)) NA else length(sg$sets$cs), susie_PP4 = NA_real_,
                      susie_PP3 = NA_real_, susie_hit_eqtl = NA_character_, susie_hit_gwas = NA_character_)
    if (!is.null(sq) && !is.null(sg) && length(sq$sets$cs) && length(sg$sets$cs)) {
      s <- coloc.susie(sq, sg)$summary
      if (!is.null(s) && nrow(s)) {
        j <- which.max(s$PP.H4.abf)
        rec$susie_PP4 <- s$PP.H4.abf[j]; rec$susie_PP3 <- s$PP.H3.abf[j]
        rec$susie_hit_eqtl <- s$hit1[j]; rec$susie_hit_gwas <- s$hit2[j]
      }
    }
    out[[length(out) + 1]] <- rec
  }
}
fwrite(rbindlist(out), file.path(dir, "results.tsv"), sep = "\t")
"""


def read_pgen(prefix, ids, samples=None):
    """Alt-allele dosages (variants x samples) for the given variant IDs from a plink2 fileset."""
    pv = pd.read_csv(f"{prefix}.pvar", sep="\t", comment="#", header=None, usecols=[2], names=["id"])
    idx = dict(zip(pv.id, range(len(pv))))
    rd = pgenlib.PgenReader(f"{prefix}.pgen".encode())
    n = rd.get_raw_sample_ct()
    out = {}
    for v in ids:
        if v in idx:
            d = np.empty(n, dtype=np.float64)
            rd.read_dosages(idx[v], d)
            d[d < 0] = np.nan
            out[v] = d if samples is None else d[samples]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", default="mage")
    ap.add_argument("--pair", nargs=2, default=["chm13v2.0", "grch38-ffx"])
    ap.add_argument("--min-pp4", type=float, default=0.5)
    args = ap.parse_args()
    a, b = args.pair
    res = ROOT / "results/06-compare/gwas"
    work = Path(os.environ["TMPDIR"]) / "coloc-susie"
    work.mkdir(parents=True, exist_ok=True)

    t = pd.read_csv(res / f"coloc-{args.cohort}-{a}-vs-{b}.tsv", sep="\t", dtype={"gene_ncbi": str})
    t = t[t[[f"PP4_{a}", f"PP4_{b}"]].max(axis=1) >= args.min_pp4]
    print(f"{len(t)} gene x trait loci with PP4 >= {args.min_pp4} on either reference", flush=True)

    m = pd.concat(pd.read_csv(f, sep="\t", usecols=["side", "category", "t2t_id", "grch38_id", "allele_match"],
                              low_memory=False)
                  for f in sorted((ROOT / "results/06-compare/variant-matching").glob("chr*.tsv.gz")))
    m = m[(m.side == "t2t") & (m.category == "matched")]
    t2g = dict(zip(m.t2t_id, m.grch38_id))
    t2flip = dict(zip(m.t2t_id, np.where(m.allele_match.eq("swapped"), -1.0, 1.0)))
    del m

    gdir = ROOT / f"results/05-eqtl/{args.cohort}/genotypes"
    n_cohort = len(pd.read_csv(gdir / "grch38.all.psam", sep="\t"))
    ped = pd.read_csv(PANEL / "20130606_g1k_3202_samples_ped_population.txt", sep=r"\s+")
    eur = ped[(ped.Superpopulation == "EUR") & (ped.FatherID == "0") & (ped.MotherID == "0")].SampleID
    eur_file = work / "eur.txt"
    eur.to_csv(eur_file, index=False, header=False)

    loci = []
    for r in t.itertuples():
        locus = f"{r.symbol}__{r.gcst}"
        d = work / locus
        d.mkdir(exist_ok=True)
        row = t.loc[r.Index]
        tops = [row[f"top_{x}"] for x in (a, b) if isinstance(row[f"top_{x}"], str)]
        if not tops:
            continue
        chrom = tops[0].split(":")[0]
        pos = [int(x.split(":")[1]) for x in tops]
        lo, hi = min(pos) - FLANK, max(pos) + FLANK

        # eQTL statistics on both references, as GRCh38 IDs
        eq = {}
        for ref in (a, b):
            gt = pd.read_csv(ROOT / f"data/star-references/{ref}/genes.tsv", sep="\t", dtype={"gene_ncbi": str})
            pids = list(gt.loc[gt.gene_ncbi == r.gene_ncbi, "gene_id"])
            parts = []
            for f in sorted(glob.glob(str(ROOT / f"results/05-eqtl/{args.cohort}/{ref}/star/matched/nominal/{chrom}.*.parquet"))):
                parts.append(pd.read_parquet(f, columns=["phenotype_id", "variant_id", "slope", "slope_se"],
                                             filters=[("phenotype_id", "in", pids)]))
            x = pd.concat(parts)
            if ref == "chm13v2.0":
                x["slope"] = x.slope * x.variant_id.map(t2flip)
                x["variant_id"] = x.variant_id.map(t2g)
                x = x.dropna(subset=["variant_id", "slope"])
            eq[ref] = x.set_index("variant_id")[["slope", "slope_se"]].rename(
                columns={"slope": f"slope_{ref}", "slope_se": f"se_{ref}"})
        v = eq[a].join(eq[b], how="inner")
        p = v.index.str.split(":")
        v["pos"] = [int(x[1]) for x in p]
        v = v[np.array([len(x[2]) == 1 and len(x[3]) == 1 for x in p]) & (v.pos >= lo).values & (v.pos <= hi).values]

        # GWAS, aligned to the variant ID's alt allele
        f = glob.glob(str(ROOT / f"data/external/gwas-sumstats/{r.gcst}/*.h.tsv.gz"))[0]
        gw = coloc.read_gwas(f, {chrom: [(lo, hi)]})
        fwd = gw.chrom + ":" + gw.pos.astype(str) + ":" + gw.a_other + ":" + gw.a_effect
        rev = gw.chrom + ":" + gw.pos.astype(str) + ":" + gw.a_effect + ":" + gw.a_other
        gb = pd.concat([pd.DataFrame({"snp": fwd, "g_beta": gw.beta, "g_se": gw.se}),
                        pd.DataFrame({"snp": rev, "g_beta": -gw.beta, "g_se": gw.se})]
                       ).drop_duplicates("snp").set_index("snp")
        v = v.join(gb, how="inner")
        meta = yaml.safe_load(open(f + "-meta.yaml"))
        gwas_n = sum(s.get("sample_size", 0) for s in meta.get("samples", []))

        # genotypes: cohort (GRCh38) and 1KGP EUR
        qd = read_pgen(str(gdir / "grch38.all"), list(v.index))
        vcf = PANEL / f"1kGP_high_coverage_Illumina.{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
        subprocess.run(["plink2", "--vcf", str(vcf), "--keep", str(eur_file), "--chr", chrom, "--from-bp", str(lo),
                        "--to-bp", str(hi), "--snps-only", "just-acgt", "--max-alleles", "2", "--set-all-var-ids",
                        "@:#:$r:$a", "--new-id-max-allele-len", "100", "--output-chr", "chrM", "--make-pgen", "--out", str(d / "eur"),
                        "--threads", "2", "--memory", "8000"], check=True, capture_output=True)
        ed = read_pgen(str(d / "eur"), list(v.index))
        keep = [s for s in v.index if s in qd and s in ed]
        Q = np.array([qd[s] for s in keep]); E = np.array([ed[s] for s in keep])
        mq = np.nanmean(Q, axis=1) / 2; me = np.nanmean(E, axis=1) / 2
        ok = (np.minimum(mq, 1 - mq) >= 0.01) & (np.minimum(me, 1 - me) >= 0.01) & ~np.isnan(Q).all(1)
        keep = [s for s, o in zip(keep, ok) if o]
        Q, E = Q[ok], E[ok]
        Q = np.where(np.isnan(Q), np.nanmean(Q, axis=1, keepdims=True), Q)
        E = np.where(np.isnan(E), np.nanmean(E, axis=1, keepdims=True), E)
        v = v.loc[keep].copy()
        v["maf_cohort"] = np.minimum(mq[ok], 1 - mq[ok]); v["maf_eur"] = np.minimum(me[ok], 1 - me[ok])
        v.index.name = "snp"
        v.reset_index().to_csv(d / "variants.tsv", sep="\t", index=False)
        pd.DataFrame(np.corrcoef(Q)).to_csv(d / "ld-cohort.tsv", sep="\t", index=False, header=False, float_format="%.5g")
        pd.DataFrame(np.corrcoef(E)).to_csv(d / "ld-eur.tsv", sep="\t", index=False, header=False, float_format="%.5g")
        loci.append({"locus": locus, "gene_ncbi": r.gene_ncbi, "symbol": r.symbol, "trait": r.trait, "gcst": r.gcst,
                     "refs": f"{a},{b}", "eqtl_n": n_cohort, "gwas_n": gwas_n, "nsnps": len(v),
                     "abf_call": r.call})
        print(f"{locus}: {len(v)} variants", flush=True)

    loci = pd.DataFrame(loci)
    loci.to_csv(work / "loci.tsv", sep="\t", index=False)
    (work / "run.R").write_text(R_SCRIPT)
    subprocess.run(["Rscript", str(work / "run.R"), str(work)], check=True)
    rs = pd.read_csv(work / "results.tsv", sep="\t")
    w = rs.pivot(index="locus", columns="ref",
                 values=["nsnps", "abf_PP4", "abf_PP3", "cs_eqtl", "cs_gwas", "susie_PP4", "susie_PP3",
                         "susie_hit_eqtl", "susie_hit_gwas"])
    w.columns = [f"{c}_{r}" for c, r in w.columns]
    out = loci.merge(w, left_on="locus", right_index=True)
    out.to_csv(res / f"coloc-susie-{args.cohort}-{a}-vs-{b}.tsv", sep="\t", index=False, float_format="%.4g")
    pd.set_option("display.width", 250)
    print(out[["symbol", "trait", "abf_call", "nsnps", f"abf_PP4_{a}", f"abf_PP4_{b}", f"cs_eqtl_{a}", f"cs_eqtl_{b}",
               f"cs_gwas_{a}", f"susie_PP4_{a}", f"susie_PP4_{b}"]].to_string(index=False))


if __name__ == "__main__":
    main()
