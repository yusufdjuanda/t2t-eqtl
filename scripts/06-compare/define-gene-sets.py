"""Step 6: the changed-gene sets used from 2026-10-29 on (decision 2026-10-29, docs 13B.4).

    pixi run python scripts/06-compare/define-gene-sets.py [--cohort mage] [--replication geuvadis]

Per reference pair, from paired-<variants>/robust-<A>-vs-<B>.tsv (paired-test.py), unflagged genes only:
  extent  robust_paired: Bonferroni-significant paired test under STAR and salmon, same sign of delta
  main    extent and |delta| >= 0.1 SD under both STAR and salmon
  core    main and also robust with the independent z (robust_indep; the original 57 for T2T vs patched GRCh38)
tier = the narrowest set the gene is in. Replication calls (same_sign, nominal, bonferroni) are added from
<replication>/paired-<variants>-at-<cohort>/replication-<A>-vs-<B>.tsv where present.
Writes results/06-compare/eqtl-comparison/<cohort>/gene-sets/changed-<A>-vs-<B>.tsv and prints the counts.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PAIRS = [("chm13v2.0", "grch38"), ("chm13v2.0", "grch38-ffx"), ("grch38-ffx", "grch38")]
FLOOR = 0.1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", default="mage")
    ap.add_argument("--replication", default="geuvadis")
    ap.add_argument("--variants", default="matched")
    args = ap.parse_args()
    base = ROOT / "results/06-compare/eqtl-comparison"
    out = base / f"{args.cohort}/gene-sets"
    out.mkdir(parents=True, exist_ok=True)
    for a, b in PAIRS:
        t = pd.read_csv(base / f"{args.cohort}/paired-{args.variants}/robust-{a}-vs-{b}.tsv", sep="\t",
                        dtype={"gene_ncbi": str})
        t = t[t.robust_paired & ~t.artefact_flag.astype(bool)].copy()
        main_ = np.minimum(t.delta.abs(), t.delta_salmon.abs()) >= FLOOR
        core = main_ & t.robust_indep
        t["tier"] = np.where(core, "core", np.where(main_, "main", "extent"))
        rep = base / f"{args.replication}/paired-{args.variants}-at-{args.cohort}/replication-{a}-vs-{b}.tsv"
        if rep.exists():
            r = pd.read_csv(rep, sep="\t", dtype={"gene_ncbi": str}).set_index("gene_ncbi")
            t = t.join(r[["same_sign", "nominal", "bonferroni"]].add_prefix(f"{args.replication}_"), on="gene_ncbi")
        cols = ["gene_ncbi", "symbol", "biotype", "region_class", "tier", "var_a", "var_b", "geno_r", "slope_a",
                "slope_b", "delta", "z_paired", "delta_salmon", "z_paired_salmon", "robust_indep"]
        cols += [c for c in t.columns if c.startswith(f"{args.replication}_")]
        t = t[cols].sort_values(["tier", "z_paired"], key=lambda s: s.map({"core": 0, "main": 1, "extent": 2})
                                if s.name == "tier" else -s.abs())
        t.to_csv(out / f"changed-{a}-vs-{b}.tsv", sep="\t", index=False, float_format="%.4g")
        n = t.tier.value_counts()
        print(f"{a} vs {b}: core {n.get('core', 0)}, main (incl. core) {n.get('core', 0) + n.get('main', 0)}, "
              f"extent (all) {len(t)}")


if __name__ == "__main__":
    main()
