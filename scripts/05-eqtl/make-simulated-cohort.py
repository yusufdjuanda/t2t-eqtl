"""Step 5 test: build a simulated cohort 'simtest' with planted eQTLs, to test the eQTL scripts end to end.

    pixi run python scripts/05-eqtl/make-simulated-cohort.py

Uses real T2T genotypes (chr22, results/01-genotypes/chm13v2.0/expression/) of 300 MAGE individuals and 300 real chm13v2.0 genes on chr22.
Expression is simulated: for 60 genes, log-expression depends on one real common variant (MAF >= 0.1, within 200 kb
of the TSS) with effect 0.3-0.8 SD per allele; 240 genes have no genetic effect. Counts ~ negative binomial.
Writes results/04-expression/simtest/matrices/chm13v2.0/{star-counts,salmon-counts,salmon-tpm}.tsv.gz,
config/samplesheets/simtest.csv (sample, sex) and results/05-eqtl/simtest/planted.tsv (gene, variant, beta).
Run the pipeline with chromosome chr22 only; a working pipeline finds most planted genes and few others.
"""
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
rng = np.random.default_rng(2026)
BCF = ROOT / "results/01-genotypes/chm13v2.0/expression/chr22.bcf"

sheet = pd.read_csv(ROOT / "config/samplesheets/mage.csv").drop_duplicates("sample")
inds = sheet["sample"].sample(300, random_state=1).tolist()
genes = pd.read_csv(ROOT / "data/star-references/chm13v2.0/genes.tsv", sep="\t", dtype={"gene_ncbi": str})
genes = genes[(genes.chrom == "NC_060946.1") & (genes.biotype == "protein_coding")].drop_duplicates("gene_id")
genes = genes.sample(300, random_state=2).reset_index(drop=True)
genes["tss"] = np.where(genes.strand == "+", genes.start, genes.end)

# common variants near the first 60 genes' TSS
q = subprocess.run(["bcftools", "query", "-s", ",".join(inds), "-i", "INFO/MAF>=0.1 && TYPE='snp'",
                    "-f", "%POS\t%REF\t%ALT[\t%GT]\n", str(BCF)], capture_output=True, text=True, check=True).stdout
rows = [l.split("\t") for l in q.strip().split("\n")]
pos = np.array([int(r[0]) for r in rows])
planted = []
dos = {}
for i in range(60):
    t = genes.tss[i]
    idx = np.where(np.abs(pos - t) <= 200_000)[0]
    if len(idx) == 0:
        continue
    j = rng.choice(idx)
    r = rows[j]
    gt = np.array([sum(int(a) for a in g.replace("|", "/").split("/") if a in "01") for g in r[3:]], dtype=float)
    vid = f"chr22:{r[0]}:{r[1]}:{r[2]}"
    beta = rng.uniform(0.3, 0.8) * rng.choice([-1, 1])
    planted.append((genes.gene_id[i], vid, beta)); dos[genes.gene_id[i]] = (gt - gt.mean()) / gt.std()

n, G = len(inds), len(genes)
base = rng.normal(5, 1.5, G)  # log mean count
anc = rng.normal(0, 1, n)     # hidden factor (like ancestry / batch)
load = rng.normal(0, 0.2, G)
logmu = base[:, None] + load[:, None] * anc[None, :] + rng.normal(0, 0.5, (G, n))
for g, v, b in planted:
    k = genes.index[genes.gene_id == g][0]
    logmu[k] += b * 0.5 * dos[g]
mu = np.exp(logmu) * rng.uniform(0.7, 1.3, n)[None, :]  # library size
r = 20
counts = rng.negative_binomial(r, r / (r + mu))
cdf = pd.DataFrame(counts, index=genes.gene_id, columns=inds)
cdf.index.name = "gene_id"
tpm = cdf / cdf.sum() * 1e6

out = ROOT / "results/04-expression/simtest/matrices/chm13v2.0"
out.mkdir(parents=True, exist_ok=True)
for name, m in [("star-counts", cdf), ("salmon-counts", cdf), ("salmon-tpm", tpm)]:
    m.round(3).to_csv(out / f"{name}.tsv.gz", sep="\t", compression="gzip")
s = sheet.set_index("sample").loc[inds].reset_index()[["sample", "sex"]]
s.to_csv(ROOT / "config/samplesheets/simtest.csv", index=False)
p = ROOT / "results/05-eqtl/simtest"; p.mkdir(parents=True, exist_ok=True)
pd.DataFrame(planted, columns=["gene_id", "variant_id", "beta"]).to_csv(p / "planted.tsv", sep="\t", index=False)
print(f"simtest: {n} individuals, {G} genes, {len(planted)} planted eQTLs")
