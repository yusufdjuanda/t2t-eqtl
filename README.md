# t2t-eqtl

**Which cis-eQTLs change when the reference genome changes?** RNA-seq of lymphoblastoid cell lines (MAGE: 731
individuals, discovery; Geuvadis: 451, replication) and 1000 Genomes genotypes, processed from the same raw reads on
three references: T2T-CHM13v2.0, GRCh38, and GRCh38 patched by FixItFelix (Behera et al. 2023). The central question is
what T2T changes *beyond* the GRCh38 errors FixItFelix already fixed.

**Status (2026-11-17).** Expression, eQTL mapping (chr1–22, STAR and salmon) and the comparison are done on MAGE and
replicated in Geuvadis. 117 genes change their eQTL between T2T and patched GRCh38 (paired test, both quantifiers),
none inside FixItFelix's regions and 74 of them in segmental duplications; 100 of 110 testable changes replicate.
Allele-specific expression, to decide which reference is right, is the next step.

| Read this | For |
|---|---|
| [`docs/t2t-eqtl.md`](docs/t2t-eqtl.md) | everything: aim, hypotheses, design, data sources, rules, how the code fits together (11.4), reproducibility set-up (11.5), step log, results |
| [`docs/data-journey/data-journey.html`](docs/data-journey/data-journey.html) | a plain-language walk through the pipeline, following two people and two genes with real numbers |
| [`config/data-sources.tsv`](config/data-sources.tsv) | every input: origin, URL or accession, version, local path, checksum |

## Requirements

- Linux x86-64 (the bioinformatics tools are not built for Windows; use WSL2 there).
- [pixi](https://pixi.sh) ≥ 0.67 (tested with 0.67.2): `curl -fsSL https://pixi.sh/install.sh | sh`
- git, and for stages 01–06 a Slurm cluster. Expression alone takes ~4,700 CPU hours for MAGE and ~2,500 for
  Geuvadis (three references); each cohort's reads (~2.5 TB) pass through scratch one batch (~160 GB) at a time.
  Inputs are ~60 GB, STAR indexes ~30 GB per reference, and `results/` ~540 GB (mostly nominal eQTL tables).
- Everything else (Python 3.11, R 4.5, tensorQTL, STAR, salmon, samtools, …) comes from `pixi.lock`.

## Quick start

```bash
git clone https://github.com/yusufdjuanda/t2t-eqtl.git && cd t2t-eqtl
pixi install --locked                 # default environment (Python, R, plink2, tensorQTL); `-e bio` for STAR & co.
cp config/site.env.example config/site.env   # edit: storage folders, scratch, Slurm partitions
pixi run setup                        # data/, results/, logs/, work/ as links (or plain folders)
pixi run fonts                        # optional: the figure font, TeX Gyre Heros
pixi run -e lint lint                 # static checks (also run by CI)
```

`config/site.env` is the only file that knows about your machine. Scripts never contain absolute paths: they find the
project root and read `config/site.env` through `scripts/lib/common.sh`.

## Running the analysis

Stages share one number across `scripts/`, `results/` and `logs/`. Slurm jobs are submitted from the project root with
`scripts/submit.sh` (or `pixi run submit`), which adds the site's partition and account:

```bash
scripts/submit.sh scripts/05-eqtl/run-tensorqtl.sbatch mage chm13v2.0
```

| Stage | What | Entry point |
|---|---|---|
| 00 | set-up: storage links, font, checksums, lint | `pixi run setup`, `fonts`, `checksums`, `verify-data`; `-e lint lint` |
| 01 | genotypes: download both 1000 Genomes panels, subset to the 1,089 expression individuals | `scripts/01-genotypes/make-sample-lists.sh`, `*/fetch-panel.*`, `*/subset-genotypes.sbatch` |
| 02 | references: STAR references and indexes for T2T, GRCh38, patched GRCh38 | `prepare-star-references.sbatch`, `prepare-ffx-reference.sbatch`, `build-star-index.sbatch` |
| 03 | five-sample pilot (go/no-go) | `scripts/03-rnaseq-pilot/` |
| 04 | expression for whole cohorts (Nextflow: SRA → STAR + salmon on each reference) | `bash scripts/04-expression/submit-cohort.sh <mage\|geuvadis>` |
| 05 | cis-eQTL mapping (tensorQTL), whole chain with dependencies | `bash scripts/05-eqtl/submit-eqtl-cohort.sh <cohort>` |
| 06 | comparison between references (variant matching, paired test, gene sets, replication, GWAS, ASE data) | `scripts/06-compare/` (order in `docs/t2t-eqtl.md` 11.4) |
| 07 | figures, tables, diagrams, workbook | `pixi run reproduce` |

Every script states its inputs, outputs and usage in its header. Downloads retry and resume, because some networks
reset HTTPS connections (see `scripts/01-genotypes/chm13v2.0/fetch-panel.sh`).

### Reproducing the figures only

The figures need only the stage 05–06 result tables, not the raw data:

```bash
pixi run reproduce                    # figure tables → figures → diagrams → workbook
pixi run -e docs data-journey         # re-render the data-journey page
```

`results/07-presentation/data/` (~5 MB) holds every table the figures are drawn from.

## Input data

Nothing is redistributed here; each source keeps its own terms. `config/data-sources.tsv` lists every input, and
`config/checksums/inputs.sha256` the sha256 of every downloaded file (check with `pixi run verify-data`). The main
inputs: 1000 Genomes phased genotypes on T2T (Zenodo 10.5281/zenodo.7612953) and on GRCh38 (NYGC 20220422); MAGE
(SRA PRJNA851328) and Geuvadis (ENA PRJEB3366) RNA-seq; NCBI RefSeq genomes and RS_2025_08 annotations; FixItFelix's
patched GRCh38 (Zenodo 7532965); GWAS Catalog summary statistics.

On the original cluster the references come from the sibling project `t2t-atlas` (`T2T_ATLAS_DIR` in
`config/site.env`); elsewhere, download them into `data/references/` as `config/data-sources.tsv` describes.

## Software environments (`pixi.toml`)

| Environment | Contents | Used by |
|---|---|---|
| `default` | Python 3.11 (pandas, scipy, PyTorch CPU, pgenlib), R (coloc, susieR, ggplot2, gt), plink2, UCSC liftOver, Graphviz, pandoc | all Python and R scripts (`pixi run python …`) |
| `bio` | STAR 2.7.11b, salmon 1.11.4, samtools 1.23.1, bcftools 1.23.1, bedtools 2.31.1, sra-tools 3.4.1, minimap2 2.31, Nextflow 25.10.4 | shell and Slurm scripts (`t2t_env bio`), Nextflow tasks |
| `longread` | Clair3 2.0.3, WhatsHap 2.8 | ONT variant calling for allele-specific expression |
| `docs` | Quarto 1.10.19 | `pixi run -e docs data-journey` |
| `lint` | Python, ShellCheck | `pixi run -e lint lint`, CI |

tensorQTL 1.0.10 and qtl 0.1.10 are vendored from their PyPI wheels in `envs/vendor/` (sha256-checked), because pixi
could not reach its PyPI mapping service from the cluster; `pixi.toml` puts them on `PYTHONPATH`.

## Conventions worth knowing

- **TMPDIR.** Every job uses a private folder under `$T2T_TMP` (`t2t_tmpdir`), never `/tmp`.
- **Chromosome names.** The panels use `chr1`; the NCBI T2T reference uses `NC_…`. Convert with
  `config/chrom-names-chm13v2.0.tsv`; never assume they match.
- **Alleles.** REF/ALT differ between T2T and GRCh38 at ~20% of common SNPs. Compared slopes are always per T2T ALT
  allele, with signs aligned by the genotype correlation (`docs/data-journey`, "Swapped alleles, step by step").
- **Resource ledger.** Every job is appended to `data/metadata/resource-usage.tsv`.

## Licence and citation

Code: MIT ([`LICENSE`](LICENSE)). Documents and figures: CC BY 4.0 ([`LICENSE-docs.md`](LICENSE-docs.md)).
Please cite as described in [`CITATION.cff`](CITATION.cff).
