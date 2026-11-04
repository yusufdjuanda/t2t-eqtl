"""Step 6: genes changed between references under both quantifications (STAR and salmon), the project's rule.

    pixi run python scripts/06-compare/combine-methods.py --cohort mage [--variants matched]

Inputs: results/06-compare/eqtl-comparison/<cohort>/{star,salmon}-<variants>/effects-<A>-vs-<B>.tsv
(compare-effects.py). A gene counts as changed for a pair when it is "changed" under both methods with the same
sign of z_diff. Writes results/06-compare/eqtl-comparison/<cohort>/robust-<variants>/{changed-<A>-vs-<B>.tsv,
summary.tsv}.

Each method picks its own lead variant pair, and in about half of the changed genes the two slopes have opposite signs,
so "stronger_on" (larger |slope|) does not describe those genes. Added columns (2026-10-27):
  pattern          "stronger on <ref>" when the STAR slopes share a sign, "opposite sign" otherwise
  salmon_z_star_var, salmon_p_star_var
                   salmon z_diff and p at STAR's own variant pair (same allele alignment sign_b = sign(geno_r)),
                   a check that both methods agree at one variant pair; not part of the selection rule
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[2]
PAIRS = [("chm13v2.0", "grch38"), ("chm13v2.0", "grch38-ffx"), ("grch38-ffx", "grch38")]


def salmon_at(cohort, variants, ref, gene_ids, var):
    """Salmon nominal slope and se for each gene (index of var) at the variant var; NaN where not tested."""
    w = pd.DataFrame({"phenotype_id": gene_ids.values, "variant_id": var.values}, index=var.index).dropna()
    w["chrom"] = w.variant_id.str.split(":").str[0]
    d = ROOT / f"results/05-eqtl/{cohort}/{ref}/salmon/{variants}/nominal"
    parts = [pd.read_parquet(d / f"{c}.cis_qtl_pairs.{c}.parquet", columns=["phenotype_id", "variant_id", "slope", "slope_se"],
                             filters=[("phenotype_id", "in", x.phenotype_id.unique().tolist())])
             for c, x in w.groupby("chrom")]
    n = pd.concat(parts).drop_duplicates(["phenotype_id", "variant_id"])
    m = w.rename_axis("gene_ncbi").reset_index().merge(n, on=["phenotype_id", "variant_id"], how="left")
    return m.set_index("gene_ncbi").reindex(var.index)[["slope", "slope_se"]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--variants", default="matched")
    args = ap.parse_args()
    base = ROOT / f"results/06-compare/eqtl-comparison/{args.cohort}"
    out = base / f"robust-{args.variants}"
    out.mkdir(parents=True, exist_ok=True)
    gid = {r: pd.read_csv(ROOT / f"data/star-references/{r}/genes.tsv", sep="\t", dtype={"gene_ncbi": str})
           .drop_duplicates("gene_ncbi").set_index("gene_ncbi").gene_id for r in ("chm13v2.0", "grch38", "grch38-ffx")}
    rows = []
    for a, b in PAIRS:
        e = {m: pd.read_csv(base / f"{m}-{args.variants}/effects-{a}-vs-{b}.tsv", sep="\t", dtype={"gene_ncbi": str})
             .set_index("gene_ncbi") for m in ("star", "salmon")}
        s, q = e["star"], e["salmon"]
        both = s.join(q[["z_diff", "changed", "stronger_on"]], rsuffix="_salmon", how="inner")
        star_changed = s[s.changed]
        robust = both[both.changed & both.changed_salmon.astype(bool) & (np.sign(both.z_diff) == np.sign(both.z_diff_salmon))]
        robust = robust.sort_values("p_diff")
        robust["pattern"] = np.where(np.sign(robust.slope_a) == np.sign(robust.slope_b),
                                     "stronger on " + robust.stronger_on, "opposite sign")
        sa = salmon_at(args.cohort, args.variants, a, gid[a].reindex(robust.index), robust.var_a)
        sb = salmon_at(args.cohort, args.variants, b, gid[b].reindex(robust.index), robust.var_b)
        robust["salmon_z_star_var"] = (sa.slope - sb.slope * np.sign(robust.geno_r)) / np.sqrt(sa.slope_se ** 2 + sb.slope_se ** 2)
        robust["salmon_p_star_var"] = 2 * norm.sf(robust.salmon_z_star_var.abs())
        robust.to_csv(out / f"changed-{a}-vs-{b}.tsv", sep="\t", float_format="%.4g")
        for so in (a, b):
            r = robust[robust.stronger_on == so]
            row = {"pair": f"{a} vs {b}", "stronger_on": so,
                   "changed_star": int((star_changed.stronger_on == so).sum()),
                   "changed_salmon": int((q.changed & (q.stronger_on == so)).sum()),
                   "changed_both": len(r), "both_no_artefact_flag": int((~r.artefact_flag).sum()),
                   "both_no_flag_opposite_sign": int((~r.artefact_flag & (r.pattern == "opposite sign")).sum())}
            row.update({f"clean_{k}": int(v) for k, v in r[~r.artefact_flag].region_class.value_counts().items()})
            rows.append(row)
        c = robust[~robust.artefact_flag]
        ok = (np.sign(c.salmon_z_star_var) == np.sign(c.z_diff)) & (c.salmon_p_star_var < 0.05 / max(len(c), 1))
        print(f"{a} vs {b}: STAR {len(star_changed)}, salmon {int(q.changed.sum())}, both {len(robust)} "
              f"({len(c)} without artefact flag, {int((c.pattern == 'opposite sign').sum())} of them opposite sign); "
              f"STAR-changed not testable with salmon: {len(star_changed.index.difference(q.index))}; "
              f"salmon at STAR's variant pair: {int(ok.sum())} of {len(c)} same sign and p < 0.05/{len(c)}")
    smry = pd.DataFrame(rows).fillna(0)
    smry.to_csv(out / "summary.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print("\n" + smry.to_string(index=False))


if __name__ == "__main__":
    main()
