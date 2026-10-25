"""Pair one chromosome's T2T and GRCh38 variants (called by match-variants.sbatch; see its header).

    python match-variants.py <work dir with t2t.tsv, g38.tsv, lifted.bed> <output .tsv.gz>
"""
import sys

import pandas as pd

COMP = str.maketrans("ACGTN", "TGCAN")
rc = lambda s: s.translate(COMP)[::-1]

d, out = sys.argv[1], sys.argv[2]
t = pd.read_csv(f"{d}/t2t.tsv", sep="\t", header=None, names=["chrom", "pos", "ref", "alt", "maf", "syntenic"],
                dtype={"chrom": str, "ref": str, "alt": str, "syntenic": str})
g = pd.read_csv(f"{d}/g38.tsv", sep="\t", header=None, names=["chrom", "pos", "ref", "alt", "maf"],
                dtype={"chrom": str, "ref": str, "alt": str})
t["row"] = range(1, len(t) + 1)
lift = pd.read_csv(f"{d}/lifted.bed", sep="\t", header=None, names=["g_chrom", "g0", "g_pos", "row", "score", "strand"],
                   dtype={"g_chrom": str})
t = t.merge(lift[["row", "g_chrom", "g_pos", "strand"]], on="row", how="left")
t["type"] = (t.ref.str.len() == 1) & (t.alt.str.len() == 1)
t["type"] = t.type.map({True: "snp", False: "indel"})
minus = t.strand == "-"
# alleles as they would read on GRCh38's forward strand (SNPs only; indels on minus-strand chains are left unmatched)
t["a1"] = t.ref.where(~minus, t.ref.map(rc))
t["a2"] = t.alt.where(~minus, t.alt.map(rc))
t.loc[minus & (t.type == "indel"), ["a1", "a2"]] = None

g["gid"] = g.chrom + ":" + g.pos.astype(str) + ":" + g.ref + ":" + g.alt
key_same = t.g_chrom + ":" + t.g_pos.astype("Int64").astype(str) + ":" + t.a1 + ":" + t.a2
key_swap = t.g_chrom + ":" + t.g_pos.astype("Int64").astype(str) + ":" + t.a2 + ":" + t.a1
gset = dict(zip(g.gid, g.maf))
t["match"] = None
t.loc[key_same.isin(gset), "match"] = "same"
t.loc[t.match.isna() & key_swap.isin(gset), "match"] = "swapped"
t["grch38_id"] = key_same.where(t.match == "same", key_swap.where(t.match == "swapped"))

# indels: the two panels often write the same indel differently (left alignment against different reference
# sequence), so unmatched indels on a plus-strand chain are paired with the nearest unmatched GRCh38 indel within
# 20 bp that has the same length change and the same MAF within 0.01 ("indel_nearby"; true pairs come from the same
# reads and individuals, so their MAFs agree closely; exact matches agree to 6 digits)
g["dlen"] = g.alt.str.len() - g.ref.str.len()
taken = set(t.grch38_id.dropna())
gi = g[(g.dlen != 0) & ~g.gid.isin(taken)].sort_values("pos")
t["dlen"] = t.alt.str.len() - t.ref.str.len()
cand = t[t.match.isna() & (t.type == "indel") & (t.strand == "+") & t.g_pos.notna()]
if len(cand) and len(gi):
    import numpy as np
    by = {k: v for k, v in gi.groupby("dlen")}
    for dl, rows in cand.groupby("dlen"):
        if dl not in by:
            continue
        gp = by[dl].pos.to_numpy(); gid = by[dl].gid.to_numpy(); gm = by[dl].maf.to_numpy()
        for i, p, mt in zip(rows.index, rows.g_pos.astype(int), rows.maf):
            j = np.searchsorted(gp, p)
            best = None
            for k in (j - 1, j):
                if 0 <= k < len(gp) and abs(gp[k] - p) <= 20 and abs(gm[k] - mt) <= 0.01 and gid[k] not in taken:
                    if best is None or abs(gp[k] - p) < abs(gp[best] - p):
                        best = k
            if best is not None:
                taken.add(gid[best]); t.loc[i, "grch38_id"] = gid[best]; t.loc[i, "match"] = "indel_nearby"
t["maf_grch38"] = t.grch38_id.map(gset)
t["category"] = "lifted_not_called_grch38"
t.loc[t.g_pos.isna(), "category"] = "no_liftover"
t.loc[t.match.notna(), "category"] = "matched"
# a true pair is the same variant in the same individuals from the same reads, so its MAF agrees (exact SNP matches
# agree to 6 digits). Pairs differing by > 0.02 are mostly paralogous positions or different indels at one site
# (genome-wide 2026-10-17: 0.24% of SNP and 5% of indel pairs); kept as their own category, not as matched.
dmaf = (t.maf - t.maf_grch38).abs()
t.loc[(t.category == "matched") & (dmaf > 0.02), "category"] = "matched_maf_discordant"
t["t2t_id"] = t.chrom + ":" + t.pos.astype(str) + ":" + t.ref + ":" + t.alt

matched = set(t.grch38_id.dropna())
g["type"] = ((g.ref.str.len() == 1) & (g.alt.str.len() == 1)).map({True: "snp", False: "indel"})
g["category"] = g.gid.isin(matched).map({True: "matched", False: "not_called_t2t"})

a = pd.DataFrame({"side": "t2t", "t2t_id": t.t2t_id, "grch38_id": t.grch38_id, "category": t.category,
                  "allele_match": t.match, "type": t.type, "maf_t2t": t.maf, "maf_grch38": t.maf_grch38,
                  "t2t_syntenic": t.syntenic.eq("1")})
b = g[g.category == "not_called_t2t"]
b = pd.DataFrame({"side": "grch38", "t2t_id": None, "grch38_id": b.gid, "category": b.category,
                  "allele_match": None, "type": b.type, "maf_t2t": None, "maf_grch38": b.maf, "t2t_syntenic": None})
pd.concat([a, b]).to_csv(out, sep="\t", index=False, compression="gzip")
