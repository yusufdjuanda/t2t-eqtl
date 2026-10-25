"""Step 6: compare eQTL effect sizes between references on the same variant (run after compare-eqtl.py).

    pixi run python scripts/06-compare/compare-effects.py --cohort mage [--method star] [--variants matched]

Why: counting eGenes present on one reference and absent on the other misses the main artefact type. PWP2 and CBS
are eGenes on every reference, but their GRCh38 signal is far too strong (docs/t2t-eqtl.md 13B.3).
For every gene tested on both references of a pair (A, B) and an eGene on at least one, each reference's lead
variant is looked up on the other reference: the same gene and the same variant (T2T <-> GRCh38 IDs from the
variant-matching tables). The slope sign is taken from the correlation of the two variants' genotypes in the cohort
(tensorQTL slopes are per ALT allele): the variant-matching label allele_match = "swapped" is wrong for some indels
written differently by the two panels (2026-10-25: e.g. T2T C:CA vs GRCh38 CA:C with genotype r = +1), which turned
identical effects into opposite-sign "changes". Pairs whose genotypes correlate with |r| < 0.9 are not compared
(different variants matched by position). Per lead variant:
    z_diff = (slope_A - slope_B) / sqrt(se_A^2 + se_B^2)
The two estimates come from the same individuals and genotypes, so they are positively correlated and this z is
conservative. The gene keeps the lead with the larger |z_diff|; "changed" = Bonferroni p_diff < 0.05 over the genes
compared. "stronger_on" = the reference with the larger |slope| at that lead.
Inputs: results/06-compare/eqtl-comparison/<cohort>/<method>-<variants>/genes.tsv, the nominal parquet files of
each reference, results/06-compare/variant-matching/. Writes effects-<A>-vs-<B>.tsv and effects-summary.tsv there.
"""
import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[2]
SIDE = {"chm13v2.0": "t2t", "grch38": "grch38", "grch38-ffx": "grch38"}
PAIRS = [("chm13v2.0", "grch38"), ("chm13v2.0", "grch38-ffx"), ("grch38-ffx", "grch38")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--method", default="star")
    ap.add_argument("--variants", default="matched")
    args = ap.parse_args()
    cdir = ROOT / f"results/06-compare/eqtl-comparison/{args.cohort}/{args.method}-{args.variants}"
    t = pd.read_csv(cdir / "genes.tsv", sep="\t", dtype={"gene_ncbi": str}, low_memory=False).set_index("gene_ncbi")

    m = pd.concat(pd.read_csv(f, sep="\t", usecols=["side", "category", "t2t_id", "grch38_id", "allele_match"],
                              low_memory=False)
                  for f in sorted(glob.glob(str(ROOT / "results/06-compare/variant-matching/chr*.tsv.gz"))))
    m = m[(m.side == "t2t") & (m.category == "matched")]
    flip = m.allele_match.eq("swapped").map({True: -1.0, False: 1.0})
    t2g = dict(zip(m.t2t_id, zip(m.grch38_id, flip))); g2t = dict(zip(m.grch38_id, zip(m.t2t_id, flip)))

    def translate(v, a, b):
        """variant ID on reference a -> (ID on reference b, sign factor)"""
        if SIDE[a] == SIDE[b]:
            return (v, 1.0)
        return (t2g if SIDE[a] == "t2t" else g2t).get(v, (None, np.nan))

    # rows to compare: gene x lead-of-reference x pair, with the variant ID and phenotype ID on each side
    rows = []
    for a, b in PAIRS:
        g = t[t[f"tested_{a}"] & t[f"tested_{b}"] & (t[f"egene_{a}"] | t[f"egene_{b}"])]
        for lead_ref, other in [(a, b), (b, a)]:
            for gene, r in g.iterrows():
                v_other, s = translate(r[f"variant_id_{lead_ref}"], lead_ref, other)
                ids = {lead_ref: r[f"variant_id_{lead_ref}"], other: v_other}
                rows.append({"pair": f"{a}|{b}", "gene_ncbi": gene, "lead_of": lead_ref, "sign_b": s,
                             "pheno_a": r[f"phenotype_id_{a}"], "pheno_b": r[f"phenotype_id_{b}"],
                             "var_a": ids[a], "var_b": ids[b]})
    q = pd.DataFrame(rows)
    # sign_b converts B's slope to A's allele coding: the sign of the genotype correlation of the two variants in the
    # cohort (identical IDs on both GRCh38 references -> +1)
    from tensorqtl import pgen
    gdir = ROOT / f"results/05-eqtl/{args.cohort}/genotypes"
    rd = {s: pgen.PgenReader(str(gdir / f"{s}.{args.variants}")) for s in ("t2t", "grch38")}
    cross = q[q.pair.map(lambda p: SIDE[p.split("|")[0]] != SIDE[p.split("|")[1]]) & q.var_a.notna() & q.var_b.notna()]
    cross = cross[["pair", "var_a", "var_b"]].drop_duplicates()
    ga_side = cross.pair.str.split("|").str[0].map(SIDE)
    geno = {}
    for s in ("t2t", "grch38"):
        v = pd.concat([cross.var_a[ga_side == s], cross.var_b[ga_side != s]]).unique().tolist()
        v = [x for x in v if x in rd[s].variant_df.index]
        geno[s] = pd.DataFrame(rd[s].read_list(v), index=v, columns=rd[s].sample_ids)
    ids = geno["t2t"].columns.intersection(geno["grch38"].columns)
    gr = {}
    for p, va, vb, sa in zip(cross.pair, cross.var_a, cross.var_b, ga_side):
        sb = "grch38" if sa == "t2t" else "t2t"
        if va in geno[sa].index and vb in geno[sb].index:
            gr[(p, va, vb)] = np.corrcoef(geno[sa].loc[va, ids].astype(float), geno[sb].loc[vb, ids].astype(float))[0, 1]
    q["geno_r"] = [1.0 if SIDE[p.split("|")[0]] == SIDE[p.split("|")[1]] else gr.get((p, va, vb), np.nan)
                   for p, va, vb in zip(q.pair, q.var_a, q.var_b)]
    print(f"label-based sign disagrees with genotype correlation: {int((np.sign(q.geno_r) != q.sign_b).sum())} of "
          f"{int(q.geno_r.notna().sum())} lookups; |r| < 0.9: {int((q.geno_r.abs() < 0.9).sum())}", flush=True)
    q["sign_b"] = np.sign(q.geno_r)
    q.loc[q.geno_r.abs() < 0.9, "sign_b"] = np.nan

    # look up slope, se, p for every needed (phenotype, variant) on each reference
    stats = {}
    for ref in sorted({r for p in PAIRS for r in p}):
        want = pd.concat([q.loc[q.pair.str.split("|").str[0] == ref, ["pheno_a", "var_a"]].set_axis(["phenotype_id", "variant_id"], axis=1),
                          q.loc[q.pair.str.split("|").str[1] == ref, ["pheno_b", "var_b"]].set_axis(["phenotype_id", "variant_id"], axis=1)]
                         ).dropna().drop_duplicates()
        parts = []
        for f in sorted(glob.glob(str(ROOT / f"results/05-eqtl/{args.cohort}/{ref}/{args.method}/{args.variants}/nominal/*.parquet"))):
            x = pd.read_parquet(f, columns=["phenotype_id", "variant_id", "slope", "slope_se", "pval_nominal"])
            parts.append(x.merge(want, on=["phenotype_id", "variant_id"]))
        stats[ref] = pd.concat(parts).drop_duplicates(["phenotype_id", "variant_id"]).set_index(["phenotype_id", "variant_id"])
        print(f"{ref}: {len(stats[ref])} of {len(want)} lookups found", flush=True)

    summ = []
    for a, b in PAIRS:
        x = q[q.pair == f"{a}|{b}"].copy()
        sa = stats[a].reindex(pd.MultiIndex.from_arrays([x.pheno_a, x.var_a]))
        sb = stats[b].reindex(pd.MultiIndex.from_arrays([x.pheno_b, x.var_b]))
        x["slope_a"], x["se_a"], x["p_a"] = sa.slope.values, sa.slope_se.values, sa.pval_nominal.values
        x["slope_b"], x["se_b"], x["p_b"] = sb.slope.values * x.sign_b.values, sb.slope_se.values, sb.pval_nominal.values
        x["z_diff"] = (x.slope_a - x.slope_b) / np.sqrt(x.se_a ** 2 + x.se_b ** 2)
        x = x.dropna(subset=["z_diff"])
        x = x.loc[x.z_diff.abs().groupby(x.gene_ncbi).idxmax()]
        n = len(x)
        x["p_diff"] = 2 * norm.sf(x.z_diff.abs())
        x["changed"] = x.p_diff * n < 0.05
        x["stronger_on"] = np.where(x.slope_a.abs() >= x.slope_b.abs(), a, b)
        x["log10p_a"], x["log10p_b"] = -np.log10(x.p_a), -np.log10(x.p_b)
        x = x.join(t[["symbol", "biotype", "region_class", "artefact_flag", f"qval_{a}", f"qval_{b}"]], on="gene_ncbi")
        cols = ["gene_ncbi", "symbol", "biotype", "region_class", "artefact_flag", "lead_of", "var_a", "var_b",
                "geno_r", "slope_a", "se_a", "log10p_a", "slope_b", "se_b", "log10p_b", "z_diff", "p_diff", "changed",
                "stronger_on", f"qval_{a}", f"qval_{b}"]
        x = x[cols].sort_values("p_diff")
        x.to_csv(cdir / f"effects-{a}-vs-{b}.tsv", sep="\t", index=False, float_format="%.4g")
        c = x[x.changed]
        for so in (a, b):
            cc = c[c.stronger_on == so]
            row = {"pair": f"{a} vs {b}", "genes_compared": n, "stronger_on": so, "changed": len(cc),
                   "changed_no_artefact_flag": int((~cc.artefact_flag).sum())}
            row.update({f"clean_{k}": int(v) for k, v in cc[~cc.artefact_flag].region_class.value_counts().items()})
            summ.append(row)
        print(f"\n{a} vs {b}: {n} genes, {len(c)} changed; top:",
              ", ".join(f"{s}(z={z:.1f},{r})" for s, z, r in zip(c.symbol.head(15), c.z_diff.head(15), c.region_class.head(15))), flush=True)
    s = pd.DataFrame(summ).fillna(0)
    s.to_csv(cdir / "effects-summary.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print("\n" + s.to_string(index=False))


if __name__ == "__main__":
    main()
