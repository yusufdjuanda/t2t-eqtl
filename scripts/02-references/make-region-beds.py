"""Stage 02 (references): write the BED of regions whose alignments are kept per sample (all other alignments are discarded).

    pixi run python scripts/02-references/make-region-beds.py

Writes data/star-references/<ref>/keep-regions.bed for chm13v2.0 and grch38, from each reference's genes.tsv:
  - each candidate gene (13.5) and, on GRCh38, its false chr21 copy, +-100 kb
  - the MHC, from the first to the last HLA gene, +-100 kb
Overlapping windows are merged. Keeping these lets us inspect the key loci later without re-aligning.
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PAD = 100_000
# Candidate loci (go/no-go 1, docs/t2t-eqtl.md 13.5 and open item 1b). Sources:
#   CBS ... KCNE1        chr21 genes next to the GRCh38 false duplication: all ten are listed in FixItFelix
#                        (Behera et al. 2023, Genome Biol, Additional file 2 Table S1); CBS, CRYAA and KCNE1 are also
#                        GIAB CMRG genes with false-duplication errors (Wagner et al. 2022, Nat Biotechnol);
#                        PWP2 and GATD3(A) are FixItFelix's artefactual eQTLs (their Table S17)
#   LOC102724560 ...     the false copies themselves: GRCh38-only RefSeq RS_2025_08 genes at chr21:5.0-6.5 Mb whose
#                        descriptions name the real gene ("cystathionine beta-synthase like" = CBS copy, etc.)
#   NOTCH2NL, SMN/NAIP, HYDIN   SD gene families resolved in T2T-CHM13; chosen from general knowledge, no single source
#   P2RY8                added after the pilot: PAR1 gene annotated on T2T chrY only (annotation artefact, 13.7)
GENES = ["CBS", "U2AF1", "CRYAA", "SIK1", "ICOSLG", "DNMT3L", "GATD3", "PWP2", "TRAPPC10", "KCNE1",
         "LOC102724560", "LOC102724594", "LOC102724652", "LOC102724428", "LOC102723996", "LOC124905051",
         "LOC102724023", "LOC102724159", "LOC102724200",
         "NOTCH2", "NOTCH2NLA", "NOTCH2NLB", "NOTCH2NLC", "NOTCH2NLR",
         "SMN1", "SMN2", "NAIP", "SERF1A", "SERF1B", "GTF2H2", "GTF2H2C", "HYDIN", "HYDIN2", "P2RY8"]

for ref in ["chm13v2.0", "grch38"]:
    g = pd.read_csv(ROOT / f"data/star-references/{ref}/genes.tsv", sep="\t")
    g = g[g.chrom.str.startswith("NC_") | ~g.chrom.str.contains("_")]  # primary chromosomes (T2T NC_…; GRCh38 no _random/_Un)
    sel = g[g.symbol.isin(GENES) | g.symbol.str.startswith("NBPF")][["chrom", "start", "end"]].copy()
    hla = g[g.symbol.str.match(r"^HLA-[A-Z]")]
    hla = hla[hla.chrom == hla.chrom.mode()[0]]
    sel = pd.concat([sel, pd.DataFrame({"chrom": [hla.chrom.iloc[0]], "start": [hla.start.min()], "end": [hla.end.max()]})])
    sel["start"] = (sel.start - 1 - PAD).clip(lower=0)
    sel["end"] = sel.end + PAD
    rows = []
    for c, df in sel.sort_values(["chrom", "start"]).groupby("chrom", sort=False):
        for s, e in zip(df.start, df.end):
            if rows and rows[-1][0] == c and s <= rows[-1][2]:
                rows[-1][2] = max(rows[-1][2], e)
            else:
                rows.append([c, s, e])
    bed = pd.DataFrame(rows, columns=["chrom", "start", "end"])
    bed.to_csv(ROOT / f"data/star-references/{ref}/keep-regions.bed", sep="\t", header=False, index=False)
    print(f"{ref}: {len(bed)} regions, {(bed.end - bed.start).sum() / 1e6:.1f} Mb; MHC {hla.chrom.iloc[0]}:"
          f"{hla.start.min()}-{hla.end.max()}")
