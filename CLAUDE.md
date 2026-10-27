# t2t-eqtl

T2T-native eQTL study: Geuvadis (and later MAGE) lymphoblastoid RNA-seq and 1KGP genotypes, processed on T2T-CHM13v2.0 and GRCh38 from the same raw reads, linked to GM12878 regulatory elements from the sibling project `~/research/t2t-atlas`. The goal is a conference submission (ASHG or ISMB 2027).

Read `docs/t2t-eqtl.md` before starting work. It is the project's only documentation file (the user asked for one file): aim, hypotheses H1–H6, design, go/no-go points, risks, execution rules (inheriting `~/research/t2t-atlas/docs/execution.md` sections 7–9), the step log, per-step details, history and open items. Add new documentation there, not in new files.

Working with the user:
- Explain concepts and results in plain language before or alongside acting. The user has asked for this explicitly.
- Keep the step log in `docs/t2t-eqtl.md` (section 12) and `data/metadata/resource-usage.tsv` current for every job.
- Do not commit unless asked.
- Do not move or rename anything in `t2t-atlas`; this project reads its references and `results/analysis/` through symlinks.

Current state (2026-10-27):
- Steps 1, 2 and 2b (pilots) are done; the user decided go (MAGE discovery, Geuvadis replication).
- Step 3: MAGE done (779 runs, 731 individuals, three references: chm13v2.0, grch38, grch38-ffx = FixItFelix patched GRCh38; matrices in `results/03-expression/mage/matrices/`). Geuvadis running (642 of 652 runs on 2026-10-27; ERR204978 has unpaired reads and was run on its pairs only, as "batch 12").
- Steps 5–6 on MAGE (chr1–22, STAR and salmon, matched variants) done: ~15,000 eGenes per reference; robust H0 set = 57 genes changed by T2T vs patched GRCh38, none in FixItFelix regions; read destination, 13 reversed eQTLs, H1 QC (~15 solid examples) done (13B.3). Step 7 presentation materials done; the manuscript is the local `.docx`.
- Next: Geuvadis replication, ctyper copy number as a covariate, chrX, allele-specific expression at the reversal loci.
- Genotypes: T2T `results/01-genotypes/chm13v2.0/expression/`, GRCh38 `results/01-genotypes/grch38/expression/` (1,089 individuals each; renamed 2026-10-27, section 16).
- `docs/t2t-eqtl.md` 11.4 explains how the scripts fit together.
- Scratch is shared with t2t-atlas and filled up on 2026-10-23 and 2026-10-26; see `docs/t2t-eqtl.md` 13A.5 and step 3c.
- Framing: FixItFelix (Behera et al. 2023) already showed artefactual GRCh38 eQTLs at the chr21 false duplication in Geuvadis. The study's central claim (H0) is what T2T changes *outside* FixItFelix's patched regions. See `docs/t2t-eqtl.md` sections 2.2, 3, 13A, 15.

Practical notes:
- HTTPS to AWS S3 and the 1000 Genomes FTP from Karakoram is reset on most TLS handshakes; use the retry/resume pattern in `scripts/01-genotypes/chm13v2.0/fetch-panel.sh`.
- Every job sets `TMPDIR` under `/scratch/$USER/t2t-eqtl/tmp/`. `/tmp` is the node's small shared root disk and was filled once on 2026-10-13.
- The panel uses `chr1`-style names; t2t-atlas `chm13v2.0-ncbi` uses `NC_…`. Map with `config/chrom-names-chm13v2.0.tsv`.
- Python: `pixi run python …` with the project's `pixi.toml`.
