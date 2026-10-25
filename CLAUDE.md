# t2t-eqtl

T2T-native eQTL study: Geuvadis (and later MAGE) lymphoblastoid RNA-seq and 1KGP genotypes, processed on T2T-CHM13v2.0 and GRCh38 from the same raw reads, linked to GM12878 regulatory elements from the sibling project `~/research/t2t-atlas`. The goal is a conference submission (ASHG or ISMB 2027).

Read `docs/t2t-eqtl.md` before starting work. It is the project's only documentation file (the user asked for one file): aim, hypotheses H1–H6, design, go/no-go points, risks, execution rules (inheriting `~/research/t2t-atlas/docs/execution.md` sections 7–9), the step log, per-step details, history and open items. Add new documentation there, not in new files.

Working with the user:
- Explain concepts and results in plain language before or alongside acting. The user has asked for this explicitly.
- Keep the step log in `docs/t2t-eqtl.md` (section 12) and `data/metadata/resource-usage.tsv` current for every job.
- Do not commit unless asked.
- Do not move or rename anything in `t2t-atlas`; this project reads its references and `results/analysis/` through symlinks.

Current state (2026-10-17):
- Steps 1, 2 and 2b (pilots) are done; the user decided go (MAGE discovery, Geuvadis replication).
- Step 3 done for MAGE: all 779 runs (731 individuals) quantified on three references (chm13v2.0, grch38, grch38-ffx = FixItFelix patched GRCh38) by 2026-10-24; matrices in `results/03-expression/mage/matrices/`. Geuvadis not yet run. Step 4 (GRCh38 genotypes), variant matching and step-5 scripts done (13B). Step 5 on MAGE (chr1–22, STAR, matched variants) done 2026-10-25 (~15,000 eGenes per reference) with STAR and salmon; robust H0 set = 57 genes changed by T2T vs patched GRCh38, none in FixItFelix regions (13B.3); next: per-gene read destination, ctyper, chrX, Geuvadis. Scratch is shared and filled up on 2026-10-23 (t2t-atlas work folder 4.4 TB); see `docs/t2t-eqtl.md` 13A.5.
- Framing: FixItFelix (Behera et al. 2023) already showed artefactual GRCh38 eQTLs at the chr21 false duplication in Geuvadis. The study's central claim (H0) is what T2T changes *outside* FixItFelix's patched regions. See `docs/t2t-eqtl.md` sections 2.2, 3, 13A, 15.

Practical notes:
- HTTPS to AWS S3 and the 1000 Genomes FTP from Karakoram is reset on most TLS handshakes; use the retry/resume pattern in `scripts/01-genotypes/fetch-1kgp-chm13-panel.sh`.
- Every job sets `TMPDIR` under `/scratch/$USER/t2t-eqtl/tmp/`. `/tmp` is the node's small shared root disk and was filled once on 2026-10-13.
- The panel uses `chr1`-style names; t2t-atlas `chm13v2.0-ncbi` uses `NC_…`. Map with `config/chrom-names-chm13v2.0.tsv`.
- Python: `pixi run python …` with the project's `pixi.toml`.
