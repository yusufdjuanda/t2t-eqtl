"""Step 6: genes changed between references under both quantifications (STAR and salmon), the project's rule.

    pixi run python scripts/06-compare/combine-methods.py --cohort mage [--variants matched]

Inputs: results/06-compare/eqtl-comparison/<cohort>/{star,salmon}-<variants>/effects-<A>-vs-<B>.tsv
(compare-effects.py). A gene counts as changed for a pair when it is "changed" under both methods with the same
sign of z_diff (the same reference shows the stronger effect). Writes
results/06-compare/eqtl-comparison/<cohort>/robust-<variants>/{changed-<A>-vs-<B>.tsv, summary.tsv}.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PAIRS = [("chm13v2.0", "grch38"), ("chm13v2.0", "grch38-ffx"), ("grch38-ffx", "grch38")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--variants", default="matched")
    args = ap.parse_args()
    base = ROOT / f"results/06-compare/eqtl-comparison/{args.cohort}"
    out = base / f"robust-{args.variants}"
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for a, b in PAIRS:
        e = {m: pd.read_csv(base / f"{m}-{args.variants}/effects-{a}-vs-{b}.tsv", sep="\t", dtype={"gene_ncbi": str})
             .set_index("gene_ncbi") for m in ("star", "salmon")}
        s, q = e["star"], e["salmon"]
        both = s.join(q[["z_diff", "changed", "stronger_on"]], rsuffix="_salmon", how="inner")
        star_changed = s[s.changed]
        robust = both[both.changed & both.changed_salmon.astype(bool) & (np.sign(both.z_diff) == np.sign(both.z_diff_salmon))]
        robust = robust.sort_values("p_diff")
        robust.to_csv(out / f"changed-{a}-vs-{b}.tsv", sep="\t", float_format="%.4g")
        for so in (a, b):
            r = robust[robust.stronger_on == so]
            row = {"pair": f"{a} vs {b}", "stronger_on": so,
                   "changed_star": int((star_changed.stronger_on == so).sum()),
                   "changed_salmon": int((q.changed & (q.stronger_on == so)).sum()),
                   "changed_both": len(r), "both_no_artefact_flag": int((~r.artefact_flag).sum())}
            row.update({f"clean_{k}": int(v) for k, v in r[~r.artefact_flag].region_class.value_counts().items()})
            rows.append(row)
        print(f"{a} vs {b}: STAR {len(star_changed)}, salmon {int(q.changed.sum())}, both {len(robust)} "
              f"({int((~robust.artefact_flag).sum())} without artefact flag); STAR-changed not testable with salmon: "
              f"{len(star_changed.index.difference(q.index))}")
    smry = pd.DataFrame(rows).fillna(0)
    smry.to_csv(out / "summary.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print("\n" + smry.to_string(index=False))


if __name__ == "__main__":
    main()
