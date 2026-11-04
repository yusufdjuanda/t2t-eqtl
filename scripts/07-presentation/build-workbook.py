"""Presentation: one Excel workbook that explains the project, its data and its first results (oral-talk format).

    pixi run python scripts/07-presentation/build-workbook.py

Inputs: results/07-presentation/data/*.tsv (prepare-figure-data.py), results/07-presentation/figures/*.png
(pipeline.dot via Graphviz; make-figures.R). Output: results/07-presentation/t2t-eqtl-project-overview.xlsx.
Text in the sheets follows docs/t2t-eqtl.md (sections 1-4, 12, 13B.3); update both together.
"""
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.drawing.image import Image
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter
from PIL import Image as PILImage

ROOT = Path(__file__).resolve().parents[2]
DD = ROOT / "results/07-presentation/data"
FD = ROOT / "results/07-presentation/figures"
OUT = ROOT / "results/07-presentation/t2t-eqtl-project-overview.xlsx"

F = "Arial"
TITLE = Font(name=F, size=16, bold=True, color="1F3864")
H2 = Font(name=F, size=12, bold=True, color="1F3864")
HEAD = Font(name=F, size=10, bold=True, color="FFFFFF")
BODY = Font(name=F, size=10)
BOLD = Font(name=F, size=10, bold=True)
NOTE = Font(name=F, size=9, italic=True, color="595959")
LINK = Font(name=F, size=10, color="0563C1", underline="single")
HEAD_FILL = PatternFill("solid", fgColor="1F3864")
BAND = PatternFill("solid", fgColor="EEF2F8")
KEY = PatternFill("solid", fgColor="FFF2CC")
THIN = Border(bottom=Side(style="thin", color="BFBFBF"))
WRAP = Alignment(wrap_text=True, vertical="top")


def table(ws, r0, c0, header, rows, widths=None, band=True):
    """Header + rows starting at (r0, c0); returns the next free row."""
    for j, h in enumerate(header):
        c = ws.cell(r0, c0 + j, h); c.font = HEAD; c.fill = HEAD_FILL; c.alignment = WRAP
    for i, row in enumerate(rows, 1):
        for j, v in enumerate(row):
            c = ws.cell(r0 + i, c0 + j, v); c.font = BODY; c.alignment = WRAP; c.border = THIN
            if band and i % 2 == 0:
                c.fill = BAND
    if widths:
        for j, w in enumerate(widths):
            ws.column_dimensions[get_column_letter(c0 + j)].width = w
    ws.freeze_panes = ws.cell(r0 + 1, c0) if r0 <= 3 else ws.freeze_panes
    return r0 + len(rows) + 2


def text(ws, r, c, v, font=BODY, merge_to=None):
    cell = ws.cell(r, c, v); cell.font = font; cell.alignment = WRAP
    if merge_to:
        ws.merge_cells(start_row=r, start_column=c, end_row=r, end_column=merge_to)
    return r + 1


def image(ws, anchor, png, width_px):
    img = Image(str(png))
    w, h = PILImage.open(png).size
    img.width, img.height = width_px, int(h * width_px / w)
    ws.add_image(img, anchor)
    return img.height


# ------------------------------------------------------------------------------------------------------------

def reversed_again(d):
    """Slopes again of opposite sign under STAR and salmon, none near zero (|slope| >= 0.05; TCAF2 and ELOB are one-sided)."""
    opp = lambda m: (np.sign(d[f"{m}_slope_a"]) != np.sign(d[f"{m}_slope_b"])) & (d[f"{m}_slope_a"].abs() >= 0.05) & (d[f"{m}_slope_b"].abs() >= 0.05)
    return opp("star") & opp("salmon")

wb = Workbook()
ws = wb.active; ws.title = "Start here"
ws.column_dimensions["A"].width = 3; ws.column_dimensions["B"].width = 30; ws.column_dimensions["C"].width = 110
r = 2
r = text(ws, r, 2, "A complete reference genome changes the eQTL map beyond known GRCh38 fixes", TITLE, 3)
r = text(ws, r, 2, f"T2T-native eQTL study in lymphoblastoid cell lines (MAGE discovery, Geuvadis replication). "
                   f"Project overview for an oral presentation (ASHG / ISMB 2027). Status as of {date.today():%d %b %Y}.", NOTE, 3)
r += 1
r = text(ws, r, 2, "What this workbook is for", H2)
r = text(ws, r, 2, "To see what the data look like, where every dataset comes from, how the pipeline runs and what the first "
                   "results say, in the order a 15-minute talk would present them. Everything here is regenerated from the "
                   "project's result files by scripts/07-presentation/ (see the Pipeline sheet).", BODY, 3)
r += 1
r = text(ws, r, 2, "The study in four sentences", H2)
for s in ["1. An eQTL is a DNA variant whose alleles go with different expression of a nearby gene; finding one depends on "
          "where RNA-seq reads are placed on the reference genome.",
          "2. GRCh38 lacks sequence, falsely duplicates some regions and collapses others; FixItFelix (Behera 2023) patched "
          "about 20 known errors and showed that some GRCh38 eQTLs (PWP2, GATD3A) are artefacts.",
          "3. We processed the same raw reads of 731 MAGE individuals on T2T-CHM13, GRCh38 and the FixItFelix-patched GRCh38, "
          "with genotypes native to each reference, and compared the eQTL effect of the same variant on the same gene.",
          "4. The patch changes 7 genes, all in its regions; T2T changes 117 more genes outside them (paired test, both "
          "quantification methods, at least 0.1 SD), mostly by moving reads between paralogs, pseudogenes and genes that exist "
          "only in T2T; they replicate in 358 independent Geuvadis individuals and change some disease colocalisations."]:
    r = text(ws, r, 2, s, BODY, 3)
r += 1
r = text(ws, r, 2, "Key numbers", H2)
key = [("Individuals (discovery, MAGE)", 731), ("Individuals (replication, Geuvadis)", 451),
       ("References compared", "T2T-CHM13v2.0, GRCh38 no-alt, FixItFelix patched GRCh38"),
       ("eGenes per reference (FDR 5%, STAR)", "15,013 / 15,259 / 15,224"),
       ("Changed genes: core (independent test) / main / extent (paired test)", "57 / 117 / 255"),
       ("Main set replicating in 358 independent Geuvadis individuals", "100 of 110 tested (86 Bonferroni)"),
       ("Genes colocalising with a GWAS trait on one reference only", "8 of 15 (3 confirmed by SuSiE-coloc)"),
       ("Genes changed by the patch (control, main definition)", 7), ("Genes changed by T2T beyond the patch (main set)", 117),
       ("... core genes with a clear read destination", "33 of 57")]
for k, v in key:
    ws.cell(r, 2, k).font = BOLD; c = ws.cell(r, 3, v); c.font = BODY; c.alignment = Alignment(horizontal="left"); r += 1
ws.cell(r, 2, "T2T effect relative to the patch").font = BOLD
c = ws.cell(r, 3, f"=ROUND(C{r-2}/C{r-3},1)&\" times as many genes as the patch fixes\""); c.font = BODY; r += 2
r = text(ws, r, 2, "Sheets", H2)
sheets = [("Talk outline", "Slide-by-slide plan of a 15-minute oral talk, with the message, the figure and speaker notes"),
          ("Data sources", "Where every dataset, reference and tool comes from: publication, accession, version, how we downloaded and checked it"),
          ("Samples", "All 1,431 sequencing runs (779 MAGE, 652 Geuvadis): individual, population, sex, depth, mapping rates"),
          ("Sample summary", "Individuals and runs per cohort and population (formulas over the Samples sheet)"),
          ("Pipeline", "Figure 1 (Graphviz) and the steps, tools, scripts and compute"),
          ("Fig2 ... Fig13", "One sheet per figure (R, ggplot2) with a caption: what it shows, how to read it, the take-home"),
          ("Key results", "The numbers behind the figures"),
          ("Gene sets", "The 255 changed genes (extent) with tier core / main / extent, paired-test changes and Geuvadis replication"),
          ("57 genes", "The 57 core genes with effect sizes and where each gene's reads go"),
          ("Glossary", "Terms used in the talk, in plain language"),
          ("Status & next", "What is done, running and next; known caveats")]
r = table(ws, r, 2, ["Sheet", "What is in it"], sheets)

# ------------------------------------------------------------------------------------------------------------
ws = wb.create_sheet("Talk outline")
ws.cell(1, 1, "Oral presentation plan (about 18 minutes + questions)").font = TITLE
slides = [
 # (slide, section: title, key message, figure, minutes, notes: why it matters / link in the story)
 (1, "Title", "A complete reference genome changes the eQTL map beyond known GRCh38 fixes", "—", 0.5,
  "Same raw reads from 731 people, three references. eQTL maps are how GWAS are read, and almost all are built on GRCh38."),
 (2, "Intro: why the reference matters", "Genotypes, read placement and gene models all depend on the reference genome", "—", 0.75,
  "GRCh38: 151 Mb unknown, false and collapsed duplications; T2T adds or corrects 238 Mb. Reads from missing sequence land "
  "elsewhere; if the loss depends on genotype, an eQTL is created or hidden."),
 (3, "Intro: gap and question", "No eQTL study has been run natively on T2T. Which eQTLs change, beyond the known patchable errors? (H0)",
  "—", 0.75, "H0: most changes lie OUTSIDE the regions FixItFelix patched; if true, patching GRCh38 is not enough."),
 (4, "Intro: has anyone done this?", "T2T used for variants, ChIP-seq, copy number, duplications, but not eQTLs; FixItFelix used T2T "
  "only to patch ~20 known GRCh38 errors", "—", 0.75,
  "Errors found by Wagner 2022 and Aganezov 2022; FixItFelix masked false copies, added missing copies from T2T, showed PWP2 "
  "and GATD3A artefacts. Its patched GRCh38 becomes our third reference."),
 (5, "Intro: the journey of the data", "1 Expression, 2 Genotypes, 3 eQTL mapping, 4 Comparison, 5 Does it matter?", "—", 0.5,
  "Design rule: same people, reads, software, annotation and covariates; only the reference changes."),
 (6, "Data: two cohorts", "731 MAGE LCLs, 26 populations (discovery); 451 Geuvadis LCLs (replication); 1000 Genomes genotypes on both references",
  "Fig2", 0.5, "MAGE: one lab, 2x150 bp. Geuvadis: several labs, 2x76 bp; 93 people in both. FixItFelix used Geuvadis, so results can be checked."),
 (7, "Design: three references", "Each run processed three times; patched vs GRCh38 checks the pipeline, T2T vs patched is the main question",
  "method-overview", 0.75, "STAR + salmon, tensorQTL with 66 covariates, plink2, liftOver variant matching; a change must appear with both quantifiers."),
 (8, "1 Expression: methods", "Download and check; STAR alignment to each reference; counts with STAR and salmon; same RefSeq genes on all", "—", 0.5,
  "STAR 2.7.11b, salmon 1.11.4; runs summed per person (731). Holds everything but the reference fixed."),
 (9, "1 Expression: mapping", "T2T halves multi-mapped reads (4.3% -> 2.1%, medians of 779 runs); the patch barely helps (4.2%)",
  "Fig3", 0.5, "STAR counts only uniquely placed reads, so T2T can count reads GRCh38 discards. H0 in miniature."),
 (10, "1 Expression: STAR and salmon", "Two ways of counting: a change counts only if both see it (T2T vs patched: 97 both, 114 STAR only, 49 salmon only)", "Fig16", 0.75,
  "STAR drops ambiguous reads; salmon shares them by likelihood. They fail differently, so requiring both removes counting artefacts; conservative (some STAR-only SD genes are probably real)."),
 (11, "1 Expression: gene level", "The patch changes 5 genes (all in its regions); T2T vs patched changes 97 of 13,563, 31% in segmental duplications",
  "Fig14", 0.75, "T2T adds no noise (99.3% within 2-fold) and shows where to look. But averages hide who loses reads."),
 (12, "2 Genotypes: methods", "Two 1000 Genomes panels; 1,089 people, frequency >= 1%; liftOver pairing; 5 genotype PCs", "—", 0.5,
  "bcftools 1.23.1, plink2 2.0.0a7.8, UCSC liftOver v482. No genotype is converted between references."),
 (13, "2 Genotypes: matched variants", "12.6 M variants paired between the T2T and GRCh38 panels; the main test uses only these",
  "Fig4", 0.5, "Same variant in the same people on both references, so any change comes from expression, not genotype calls."),
 (14, "3 eQTL mapping: methods", "TMM + inverse-normal expression; 66 covariates; tensorQTL within 1 Mb; 10,000 permutations, q < 0.05", "—", 0.5,
  "tensorQTL 1.0.10. Same recipe on every reference and counting method."),
 (15, "3 eQTL mapping: overall", "About 15,000 eGenes on every reference: the overall map is stable", "Fig5", 0.5,
  "Pipeline sanity check (agrees with MAGE paper). Totals hide which genes change."),
 (16, "3 eQTL mapping: expression is not eQTL", "Of 59 genes with a 2-fold expression change only 15 change their eQTL; 89 of 104 eQTL "
  "changes come with < 2-fold expression change", "Fig15", 1.0,
  "PKD1P4-NPIPA8: uniform loss, eQTL unchanged. CCZ1: same mean, eQTL reverses. Justifies comparing the same variant's effect."),
 (17, "4 Comparison: methods", "Same gene, same variant; paired test; a change needs all four: real (Bonferroni, |z| > 4.6), big enough (>= 0.1 SD; typical eQTL 0.30), STAR and salmon same direction, no flag", "—", 0.75,
  "Three pairs of references; region labels; checks: Geuvadis replication and copy number."),
 (18, "4 Comparison: how many change", "T2T vs GRCh38 121 (6 in patched regions); the patch alone 7, all inside; T2T vs patched 117, none inside",
  "Fig20", 1.0, "Central result (H0) with its sanity check. Robust: 76 at 0.2 SD; 100 of 110 replicate in Geuvadis."),
 (19, "4 Comparison: top 10 changed genes", "R gt table: gene, rsID, GRCh38 position, REF/ALT, effect allele and its 1000 Genomes frequency, effects, change, paired z and P, replication",
  "Fig23", 1.0, "Names and numbers early, so the next slides explain what the audience has seen (ASAH2, SMN2, CCZ1, PMS2CL)."),
 (20, "4 Comparison: segmental duplications", "50% of changed genes lie in segmental duplications vs 2.6% of all genes compared",
  "Fig21", 0.75, "Short reads cannot tell near-identical copies apart; the reference decides where they go."),
 (21, "4 Comparison: reads move (coverage)", "RNA coverage along ASAH2 and ASAH2B, 0 vs 2 copies of the eQTL allele, on both references: GRCh38 ASAH2 0.47x / ASAH2B 1.40x; T2T ASAH2 4.5x / ASAH2B 0.87x",
  "Fig25", 1.0, "Genome-browser view from the ASE target BAMs (347 people, unique reads). ASAH2 is the largest change; CCZ1, PMS2CL and SMN2 were weaker or had too few people with 0 copies."),
 (22, "4 Comparison: the eQTL reverses (ASAH2)", "ASAH2 loses what ASAH2B gains per person (r = -0.91), and the ASAH2 eQTL reverses (+0.77 vs -0.53)",
  "Fig22", 1.0, "Mechanism: genotype-dependent read placement between copies. Same in ASAH2/B, PMS2/PMS2CL, RPSA/RPSA2."),
 (23, "5 Link to GWAS", "Colocalisation with six GWAS on each reference: 6 changed genes colocalise on one reference only",
  "Fig24", 1.0, "IFITM3-AS1 with hypothyroidism on patched GRCh38 only; 16p11.2 SMG1P6 vs NPIPB13."),
 (24, "Limits and next steps", "Which reference is right (allele-specific expression running); salmon alone over-calls; unmapped reads; chrX; H5", "—", 0.75,
  "Each caveat has a check running or planned; none affects H0."),
 (25, "Conclusions", "A complete reference changes eQTLs well beyond known patches, through read placement in duplicated and new "
  "sequence, and this changes disease links", "—", 0.5,
  "Take-home: check gene families on a complete reference before trusting their eQTLs. Code and results will be released."),
 ("B1", "Backup: key numbers", "eGenes and changed genes per comparison", "—", 0, "For questions."),
 ("B2", "Backup: one gene, two people (PWP2)", "On GRCh38 the read loss depends on genotype, so the eQTL is four times too strong", "—", 0, "For questions; slide 13 makes the same point with new genes."),
 ("B3", "Backup: question and mechanisms", "H0 and H1-H4", "—", 0, "For questions."),
 ("B4", "Backup: detailed pipeline", "Download, quantify, genotypes, eQTLs, compare", "Fig1", 0, "For questions."),
 ("B5", "Backup: data sources", "Public data only; every download verified", "—", 0, "For questions."),
 ("B6", "Backup: how the tools work", "STAR, salmon, phenotypes, tensorQTL, variant pairs, paired test", "—", 0, "For questions."),
 ("B7", "Backup: replication in Geuvadis", "358 individuals not in MAGE: 100 of 110 replicate (86 Bonferroni); 10 of 12 reversals", "Fig12", 0, "For questions."),
 ("B8", "Backup: copy number", "ctyper copy number as a covariate leaves 55 of 62 covered genes changed", "Fig18", 0, "For questions."),
 ("B9", "Backup: examples", "PWP2, CCZ1 (reversed), SMN2 by genotype on two references", "Fig11", 0, "For questions."),
 ("B10", "Backup: allele-specific expression", "Which reference is right: within-person test, running", "—", 0, "For questions."),
 ("B11", "Backup: why a paired test", "Both effects come from the same people; the standard test is 9-15x too conservative",
  "Fig8", 0, "Residual correlation 0.99; covariance of correlated regressions (Zellner 1962). Bootstrap 0.95x paired se; null lambda 1.03."),
 ("B12", "Backup: how the paired test is calculated", "One gene (IFITM2): per-person expression on both references, the per-person difference by genotype, the formula and references", "Fig19", 0, "For questions; Zellner 1962, Weesie 1999, Lin & Sullivan 2009."),
 ("B13", "Backup: ASAH2 coverage, all genotypes", "Fig 25 with heterozygotes: 57 A/A, 140 A/C, 150 C/C; heterozygotes in between", "Fig25s", 0, "For questions; 347 re-aligned people."),
]
r = table(ws, 3, 1, ["Slide", "Section: title", "Key message (one sentence on the slide)", "Figure", "Minutes", "Why it matters / link in the story"],
          slides, widths=[7, 24, 60, 9, 9, 90])
ws.cell(r, 4, "Total").font = BOLD
c = ws.cell(r, 5, f"=SUM(E4:E{r-2})"); c.font = BOLD

# ------------------------------------------------------------------------------------------------------------
ws = wb.create_sheet("Data sources")
ws.cell(1, 1, "Where every dataset comes from").font = TITLE
ws.cell(2, 1, "Full details and exact download commands: docs/t2t-eqtl.md section 4. Paths are relative to the project "
              "root; data/ lives on /storage/projects/yusuf.maulana/t2t-eqtl/.").font = NOTE
src = [
 ("RNA-seq", "MAGE RNA-seq (discovery cohort)", "LCL RNA-seq of 731 1000 Genomes individuals, 26 populations, one lab; 779 runs, 2x150 bp, unstranded, NovaSeq 6000, median 28.4 M read pairs",
  "Taylor et al. 2024, Nature (MAGE)", "NCBI SRA PRJNA851328", "run report 2026-10-17", "2.51 TB reads (transient on scratch)",
  "scripts/04-expression/fetch-sra-batch.sbatch (prefetch)", "vdb-validate (NCBI checksums) + read pairs equal to ENA", "Expression on all three references"),
 ("RNA-seq", "Geuvadis RNA-seq (replication cohort)", "LCL RNA-seq of 462 individuals, 5 populations, 7 labs; 652 runs used (451 individuals with T2T genotypes), 2x76 bp, unstranded",
  "Lappalainen et al. 2013, Nature", "ENA PRJEB3366 (mirrored at NCBI SRA); sample sheet E-GEUV-1 (BioStudies)", "2026-10-17", "~2.5 TB reads",
  "scripts/04-expression/fetch-sra-batch.sbatch", "vdb-validate + read pairs equal to ENA", "Replication (done 2026-10-27)"),
 ("Metadata", "Run reports and sample sheets", "Run, individual, population, sex, lab, read counts",
  "ENA portal API; EBI BioStudies", "data/metadata/{mage,geuvadis}/", "2026-10-13 / 17", "small",
  "curl (commands in doc 4.3)", "—", "config/samplesheets/{mage,geuvadis}.csv"),
 ("Genotypes", "1000 Genomes on T2T-CHM13v2.0 (phased)", "3,202 individuals, SNVs and indels called on CHM13, phased with SHAPEIT5 v1.1",
  "Aganezov 2022 Science; Rhie 2023 Nature; Lalli et al. (panel)", "Zenodo 10.5281/zenodo.7612953; s3://human-pangenomics/.../Phased_SHAPEIT5_v1.1/", "2026-10-13", "13.0 GB BCF",
  "scripts/01-genotypes/chm13v2.0/fetch-panel.sh (curl with retry/resume)", "file sizes equal S3", "T2T-native genotypes"),
 ("Genotypes", "1000 Genomes on GRCh38 (NYGC phased panel)", "3,202 individuals, SNV/indel/SV, high-coverage GRCh38 calls",
  "Byrska-Bishop et al. 2022, Cell", "ftp.1000genomes.ebi.ac.uk .../1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/", "2026-10-17", "~30 GB",
  "scripts/01-genotypes/grch38/fetch-panel.sbatch", "md5 against 20220804_manifest.txt", "GRCh38 genotypes (both GRCh38 references)"),
 ("Reference", "T2T-CHM13v2.0 genome", "Complete human genome (24 chromosomes, chrY from HG002); chrY PARs masked, chrM and chrEBV added for RNA-seq",
  "Nurk et al. 2022, Science", "NCBI GCF_009914755.1", "2026-09-27 (cluster bundle)", "3.1 Gb",
  "t2t-atlas prepare-references.sbatch; scripts/02-references/prepare-star-references.sbatch", "CHECKSUMS.md5", "Reference 1"),
 ("Reference", "GRCh38 no-alt analysis set", "Primary assembly + unplaced contigs, rCRS chrM, chrEBV, chrY PARs masked (the set ENCODE uses)",
  "Genome Reference Consortium / NCBI", "NCBI GCA_000001405.15 seqs_for_alignment_pipelines", "2026-10-04", "3.1 Gb",
  "t2t-atlas prepare-noalt-reference.sbatch", "NCBI md5checksums.txt", "Reference 2"),
 ("Reference", "FixItFelix patched GRCh38", "GRCh38 no-alt with false duplications masked and decoys for collapsed genes",
  "Behera et al. 2023, Genome Biology", "Zenodo 7532965 (GRCh38_masked_v2_decoy_gene.fasta.gz)", "2026-10-17", "873 MB",
  "scripts/02-references/prepare-ffx-reference.sbatch", "md5 35d00c3c48d1849adb8e02086c3a1c70", "Reference 3 (separates patch effects from completeness)"),
 ("Annotation", "NCBI RefSeq gene models RS_2025_08", "Same annotation release on both assemblies: T2T 58,563 genes, GRCh38 59,257; matched by NCBI GeneID (55,334 on both)",
  "NCBI RefSeq", "GCF_009914755.1-RS_2025_08; GCF_000001405.40-RS_2025_08", "2026-09/10", "—",
  "cluster bundle; registered in t2t-atlas", "—", "Gene counts, gene matching between references"),
 ("Region tracks", "T2T sequence absent from GRCh38", "616 regions, 251.3 Mb (chm13v2-unique_to_hg38.bed)",
  "T2T consortium", "s3://human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/", "2026-10-17", "small",
  "curl with retry", "local md5 63b48dfa…", "Gene class 't2t_new_sequence'"),
 ("Region tracks", "Segmental duplications, satellites, RepeatMasker (T2T)", "SD track with identity colours; CenSat v2.1; RepeatMasker 4.1.2p1",
  "T2T consortium (Vollger 2022 Science for SDs)", "s3://human-pangenomics/T2T/CHM13/assemblies/annotation/", "2026-10-06", "small",
  "t2t-atlas", "CHECKSUMS.md5", "Gene class 'sd' (>= 98% identity)"),
 ("Region tracks", "FixItFelix regions and results", "Tables S1-S19: falsely duplicated and collapsed regions, eQTL results",
  "Behera et al. 2023, Genome Biology", "Additional file 2 (Springer static content)", "2026-10-17", "xlsx",
  "curl", "—", "Gene class 'ffx_patched'; positive control"),
 ("Region tracks", "UCSC liftOver chains (hs1 <-> hg38)", "Used only to PAIR variants called independently on each reference",
  "UCSC Genome Browser", "hgdownload.soe.ucsc.edu goldenPath hs1/hg38 liftOver", "2026-10-04", "small",
  "t2t-atlas", "CHECKSUMS.md5", "Variant matching"),
 ("Prior work", "ctyper paralog-specific copy numbers for 1KGP", "3,202 samples x 84,346 pangenome-allele types; Geuvadis expression matrix",
  "Ma & Chaisson 2025, Nature Genetics", "Zenodo 10.5281/zenodo.13683015", "2026-10-23", "tar.gz",
  "curl", "md5 equal to Zenodo", "Copy-number covariate (paired-test.py --copy-number)"),
 ("GWAS", "GWAS Catalog associations", "1.2 M associations, all studies (release 'latest')",
  "Sollis et al. 2023, Nucleic Acids Research", "ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/", "2026-10-29", "~0.3 GB",
  "curl with retry", "CHECKSUMS.md5 (local)", "GWAS lookup (gwas-lookup.py)"),
 ("GWAS", "Six GWAS summary statistics (harmonised GRCh38)", "Hypothyroidism GCST90627749 (White 2026); SLE GCST003156 (Bentham 2015); "
  "ulcerative colitis GCST004133 (de Lange 2017); multiple sclerosis GCST90558092 (Shigesi 2025); type 2 diabetes GCST006867 "
  "(Xue 2018); atrial fibrillation GCST90624412 (Yuan 2025)", "see column 'What it is'", "GWAS Catalog summary statistics FTP",
  "2026-10-29", "~4 GB", "curl with retry", "CHECKSUMS.md5 (local)", "Colocalisation (coloc.py, coloc-susie.py)"),
 ("Long reads", "1000 Genomes ONT alignments to T2T", "Per-person CRAMs (260 MAGE individuals); only target regions streamed",
  "Schloissnig et al. 2025, Nature", "ftp.1000genomes.ebi.ac.uk .../1KG_ONT_VIENNA/t2t/", "2026-10-29", "~10 Mb per person",
  "samtools view -L (ont-call.sbatch)", "retried reads; Clair3 PASS calls", "SNPs in SD genes for allele-specific expression"),
 ("Sibling project", "GM12878 ChIP-seq peaks on both references", "ENCODE GM12878 ChIP-seq reprocessed from raw reads on T2T and GRCh38",
  "ENCODE; project t2t-atlas", "data/external/t2t-atlas-analysis (read-only link)", "—", "—",
  "t2t-atlas", "—", "Planned: regulatory elements at changed eQTLs (H5)"),
 ("Software", "STAR 2.7.11b, salmon 1.11.4, sra-tools 3.4.1", "Alignment and quantification",
  "Dobin 2013; Patro 2017", "cluster modules / bioconda", "—", "—", "Nextflow module quantify_run", "—", "Step 3"),
 ("Software", "tensorQTL 1.0.10 (+ qtl 0.1.10), plink2, bcftools 1.23.1", "cis-eQTL mapping with permutations; genotype handling",
  "Taylor-Weiner et al. 2019, Genome Biology", "PyPI wheels (sha256 checked) vendored in envs/vendor/", "2026-10-17", "—",
  "pixi.toml", "sha256", "Steps 4-5"),
 ("Software", "Nextflow DSL2, pixi, R 4.5 (ggplot2), Graphviz 14", "Workflow, environments, figures", "—", "conda-forge / bioconda", "—", "—",
  "pixi.toml", "—", "Pipeline and this workbook"),
 ("Software", "R coloc 5.2.3, susieR 0.14; Clair3 2.0.3, WhatsHap 2.8, minimap2 2.31", "Colocalisation; long-read SNP calling and "
  "phasing; mappability of SNP windows", "Giambartolomei 2014; Wallace 2021; Wang 2020", "conda-forge / bioconda; cluster module", "2026-10-29",
  "—", "pixi.toml (default and 'longread' environments)", "—", "Steps 6 (H6, allele-specific expression)"),
]
r = table(ws, 4, 1, ["Type", "Resource", "What it is", "Publication / origin", "Accession / URL", "Retrieved", "Size",
                     "How obtained (script)", "How checked", "Used for"], src,
          widths=[12, 30, 52, 30, 45, 14, 14, 40, 30, 32])

# ------------------------------------------------------------------------------------------------------------
co = pd.read_csv(DD / "cohorts.tsv", sep="\t")
qc = pd.read_csv(DD / "mapping-qc.tsv", sep="\t")[["run", "unique_pct_chm13v2.0", "unique_pct_grch38", "unique_pct_grch38-ffx",
                                                   "multi_pct_chm13v2.0", "multi_pct_grch38", "multi_pct_grch38-ffx"]]
co = co.merge(qc, on="run", how="left")
co["first_run"] = (~co.duplicated(["cohort", "sample"])).astype(int)
co["cohort"] = co.cohort.map({"mage": "MAGE", "geuvadis": "Geuvadis"})
cols = ["cohort", "sample", "run", "population", "superpop", "sex", "read_pairs", "read_length", "lab", "first_run",
        "unique_pct_chm13v2.0", "unique_pct_grch38", "unique_pct_grch38-ffx", "multi_pct_chm13v2.0", "multi_pct_grch38", "multi_pct_grch38-ffx"]
ws = wb.create_sheet("Samples")
hdr = ["Cohort", "Individual", "Run", "Population", "Super-pop.", "Sex", "Read pairs", "Read length", "Lab",
       "First run of individual (1/0)", "Unique % T2T", "Unique % GRCh38", "Unique % patched", "Multi % T2T", "Multi % GRCh38", "Multi % patched"]
rows = [[None if pd.isna(v) else v for v in rec] for rec in co[cols].itertuples(index=False)]
for j, h in enumerate(hdr, 1):
    c = ws.cell(1, j, h); c.font = HEAD; c.fill = HEAD_FILL; c.alignment = WRAP
for i, row in enumerate(rows, 2):
    for j, v in enumerate(row, 1):
        ws.cell(i, j, v).font = BODY
for j, w in enumerate([10, 11, 13, 11, 10, 8, 12, 10, 10, 12, 10, 10, 10, 10, 10, 10], 1):
    ws.column_dimensions[get_column_letter(j)].width = w
for row in ws.iter_rows(min_row=2, min_col=7, max_col=7):
    for c in row:
        c.number_format = "#,##0"
ws.freeze_panes = "C2"; ws.auto_filter.ref = f"A1:P{len(rows)+1}"
n_last = len(rows) + 1
note_r = n_last + 2
ws.cell(note_r, 1, "Mapping rates (STAR) are filled for MAGE; Geuvadis rates are in results/04-expression/geuvadis/matrices/runs.tsv. Two Geuvadis runs "
                   "(NA07000, HG00237) have no population in E-GEUV-1. Source: config/samplesheets/*.csv, "
                   "results/04-expression/mage/matrices/runs.tsv.").font = NOTE

# sample summary with formulas over the Samples sheet
ws = wb.create_sheet("Sample summary")
ws.cell(1, 1, "Individuals and runs per cohort and population").font = TITLE
ws.cell(2, 1, "All numbers are formulas over the Samples sheet (individuals = runs flagged 'first run of individual').").font = NOTE
pops = co.drop_duplicates(["cohort", "population"]).dropna(subset=["population"]).sort_values(["cohort", "superpop", "population"], ascending=[False, True, True])
hdr = ["Cohort", "Super-pop.", "Population", "Individuals", "Runs"]
for j, h in enumerate(hdr, 1):
    c = ws.cell(4, j, h); c.font = HEAD; c.fill = HEAD_FILL
rng = lambda col: f"Samples!${col}$2:${col}${n_last}"
r = 5
for rec in pops.itertuples():
    ws.cell(r, 1, rec.cohort); ws.cell(r, 2, rec.superpop); ws.cell(r, 3, rec.population)
    ws.cell(r, 4, f"=SUMIFS({rng('J')},{rng('A')},A{r},{rng('D')},C{r})")
    ws.cell(r, 5, f"=COUNTIFS({rng('A')},A{r},{rng('D')},C{r})")
    for j in range(1, 6):
        ws.cell(r, j).font = BODY
    r += 1
last = r - 1
r += 1
for coh in ("MAGE", "Geuvadis"):
    ws.cell(r, 1, coh).font = BOLD; ws.cell(r, 3, "All").font = BOLD
    ws.cell(r, 4, f"=SUMIFS({rng('J')},{rng('A')},A{r})").font = BOLD
    ws.cell(r, 5, f"=COUNTIFS({rng('A')},A{r})").font = BOLD
    for c in (ws.cell(r, j) for j in range(1, 6)):
        c.fill = KEY
    r += 1
ws.cell(r, 1, "Both cohorts").font = BOLD
ws.cell(r, 4, f"=SUM(D{r-2}:D{r-1})").font = BOLD; ws.cell(r, 5, f"=SUM(E{r-2}:E{r-1})").font = BOLD
ws.cell(r + 1, 1, "93 individuals are in both cohorts (counted in each).").font = NOTE
for j, w in enumerate([12, 11, 12, 12, 10, 18], 1):
    ws.column_dimensions[get_column_letter(j)].width = w

# ------------------------------------------------------------------------------------------------------------
ws = wb.create_sheet("Pipeline")
ws.cell(1, 1, "Figure 1. Data sources and pipeline (Graphviz: scripts/07-presentation/pipeline.dot)").font = TITLE
ws.column_dimensions["A"].width = 2
h = image(ws, "B3", FD / "fig01-pipeline.png", 900)
r = 3 + int(h / 20) + 2
steps = [
 ("1", "Genotypes on T2T", "Subset the phased 1KGP T2T panel to the expression individuals", "bcftools 1.23.1", "scripts/01-genotypes/", "minutes"),
 ("2 / 2b", "Pilots", "5 Geuvadis + 5 MAGE runs on both references; go/no-go 1 passed", "STAR, salmon", "scripts/03-rnaseq-pilot/", "~2 h"),
 ("3", "Expression (Nextflow)", "Per run: SRA download and checks, STAR GeneCounts and salmon on three references, read counts per contig; "
  "merged into gene x individual matrices", "sra-tools, STAR 2.7.11b, salmon 1.11.4, Nextflow", "scripts/04-expression/ (nextflow/, submit-cohort.sh)",
  "MAGE: 6,157 CPU h (779 runs x 3 references)"),
 ("4", "Genotypes on GRCh38", "Download and subset the NYGC phased panel", "curl, bcftools", "scripts/01-genotypes/grch38/", "hours"),
 ("6a", "Variant matching", "Pair T2T and GRCh38 variants (liftOver only to pair, same alleles, MAF within 0.02); 12.6 M matched", "UCSC liftOver, pandas",
  "scripts/06-compare/match-variants.*", "~1 h"),
 ("5", "cis-eQTL mapping", "TMM + inverse-normal expression; 5 genotype PCs, sex, 60 expression PCs; map_cis (10,000 permutations) and map_nominal per chromosome",
  "tensorQTL 1.0.10, plink2", "scripts/05-eqtl/", "~44 CPU h per reference and method"),
 ("6", "Comparison", "eGene presence, same-variant effect sizes (z), STAR AND salmon rule, region classes and artefact flags, read destination",
  "pandas, scipy", "scripts/06-compare/", "minutes"),
 ("6b", "Paired test and gene sets", "Effect difference with the covariance of the two slopes from the same individuals; bootstrap "
  "calibration; gene sets core / main / extent", "numpy, pandas", "scripts/06-compare/paired-test.py, define-gene-sets.py", "~17 min per run"),
 ("6c", "Replication and copy number", "Paired test in the 358 Geuvadis individuals not in MAGE at MAGE's variant pairs; ctyper copy "
  "number as covariate", "numpy, pandas", "scripts/06-compare/replicate-changes.py, ctyper-copy-number.py", "minutes"),
 ("6d", "GWAS", "GWAS Catalog lookup with matched background; coloc on both references; SuSiE-coloc with LD", "coloc, susieR, plink2",
  "scripts/06-compare/gwas-lookup.py, coloc.py, coloc-susie.py", "38 min; 2 h"),
 ("6e", "Allele-specific expression (running)", "Feasibility; RNA realignment of target genes in 347 individuals; ONT SNP calls in 260",
  "STAR, minimap2, Clair3, WhatsHap", "scripts/06-compare/ase-*.py, ase-realign.sbatch, ont-call.sbatch", "~24 min per run; ~15 min per person"),
 ("7", "Presentation", "Figure data, R figures, Graphviz diagram, this workbook", "R 4.5 ggplot2, Graphviz, openpyxl", "scripts/07-presentation/", "minutes"),
]
ws.cell(r, 2, "Steps").font = H2
table(ws, r + 1, 2, ["Step", "Name", "What happens", "Tools", "Scripts", "Compute"], steps, widths=[8, 22, 70, 34, 38, 30])

# ------------------------------------------------------------------------------------------------------------
figs = [
 ("Fig2 Cohorts", "fig02-cohorts.png", "Figure 2. The two cohorts",
  "Top: individuals per 1000 Genomes population (colour = continental super-population). Bottom: read pairs per sequencing run.",
  "MAGE covers all 26 populations evenly (~28 each); Geuvadis is 4 European populations plus Yoruba. MAGE runs are deeper and "
  "twice as long (2x150 bp), which helps place reads in duplicated regions.",
  "MAGE is the discovery cohort because it is larger, more diverse and sequenced in one lab."),
 ("Fig3 Mapping", "fig03-mapping-qc.png", "Figure 3. Read mapping on the three references",
  "Each grey line is one MAGE run aligned to the three references; boxes summarise the 779 runs.",
  "Uniquely mapped reads rise slightly on T2T (median 84.9% vs 84.5%) and multi-mapped reads halve (2.1% vs 4.2%).",
  "The same reads are placed less ambiguously on a complete reference."),
 ("Fig4 Variants", "fig04-variant-matching.png", "Figure 4. Pairing variants between the two genotype panels",
  "Common variants (MAF >= 1%) in the 1,089 expression individuals, by whether a T2T variant has a GRCh38 partner.",
  "12.6 M variants are matched. 0.13 M pairs differ in allele frequency (paralogous positions) and are excluded. "
  "3.75 M GRCh38 variants have no T2T partner (mostly filtered by the T2T panel's QC; some are paralog artefacts).",
  "The main comparison uses matched variants only, so a change cannot come from different genotypes."),
 ("Fig5 eGenes", "fig05-egenes.png", "Figure 5. eGenes per reference",
  "Genes with a significant cis-eQTL (FDR 5%) on each reference, with STAR and with salmon counts.",
  "About 15,000 of ~16,000-17,000 tested genes on every reference, as in the MAGE paper.",
  "The overall eQTL map is stable; the reference matters for specific genes."),
 ("Fig6 Control", "fig06-positive-control.png", "Figure 6. Positive control: FixItFelix's chr21 genes",
  "Significance (-log10 q, log scale) of five genes at the chr21 false duplication on each reference.",
  "All five are far more significant on unpatched GRCh38 (PWP2 q = 1e-210) than on T2T or the patched GRCh38 (~1e-33), "
  "which agree with each other.",
  "Our pipeline reproduces the published artefacts; T2T fixes them as the patch does."),
 ("Fig7 Effects", "fig07-effect-comparison.png", "Figure 7. Same gene, same variant on T2T vs patched GRCh38",
  "Each point is a gene: the effect (slope) of its lead variant on patched GRCh38 (x) and on T2T (y). On the diagonal = no change. "
  "Colour: core (old test), the further main genes (paired test, at least 0.1 SD), extent only (smaller changes).",
  "Almost all of ~15,000 genes lie on the diagonal. 117 change by at least 0.1 SD with both quantifiers and no artefact flag "
  "(255 including smaller changes). None lies in a FixItFelix region; half are in segmental duplications.",
  "T2T changes eQTLs beyond the known patches (H0)."),
 ("Fig8 Paired test", "fig08-paired-test.png", "Figure 8. The paired test",
  "Left: z of the effect difference with the old independent-samples test (x) and the paired test (y), every compared gene; "
  "dashed: Bonferroni. Top right: bootstrap SD of the difference divided by each formula's se (1 = correct), with genomic "
  "inflation lambda. Bottom right: share replicating in 358 independent Geuvadis individuals for each definition of a change.",
  "The independent test is 9-15 times too conservative; the paired test matches the bootstrap (0.95) and has lambda 1.03 on the "
  "null pair. Replication rises with the size floor and levels off at 0.1 SD (91% of the main set).",
  "Comparing references on the same people needs a paired test; the main set is 117 genes."),
 ("Fig9 Core genes", "fig09-robust-genes.png", "Figure 9. The 57 core genes (old independent test; mechanism analyses)",
  "Left: effect difference z per gene, independent test (colour = region class). Right top: effect pattern at the same variant (stronger on "
  "one reference with the same sign; reversed = opposite signs, both p < 1e-5; opposite sign with one effect weak). "
  "Right bottom: how strongly each gene's read loss is mirrored by one partner gene.",
  "13 stronger on T2T, 13 on patched GRCh38, 13 reversed, 18 opposite sign with one weak effect. 33/57 have a clear partner "
  "gene that takes their reads.",
  "The changes are a read-placement phenomenon in gene families."),
 ("Fig10 Read swaps", "fig10-paralog-swaps.png", "Figure 10. Where the reads go",
  "Each dot is one person: change in the gene's expression (T2T minus patched GRCh38) vs change in its partner.",
  "CCZ1/CCZ1B and ASAH2/ASAH2B swap reads almost exactly; RPSA loses reads to its processed pseudogene RPSA2; "
  "PDPR on GRCh38 holds part of the reads that T2T assigns to a gene in sequence GRCh38 lacks (r = -0.82, about 8% of its change).",
  "Missing or collapsed sequence moves reads between copies, and the move differs between people."),
 ("Fig11 Examples", "fig11-eqtl-examples.png", "Figure 11. Three eQTLs on two references",
  "Normalised expression by genotype at the lead variant, for the same people, on two references.",
  "PWP2: strong genotype effect only on unpatched GRCh38 (the known artefact). CCZ1: the effect reverses direction. "
  "SMN2: a much stronger effect on T2T.",
  "The reference can create, remove or reverse an eQTL."),
 ("Fig12 Replication", "fig12-replication.png", "Figure 12. Replication in Geuvadis at the same variants",
  "Left: paired z of the effect difference (T2T minus patched GRCh38) for each main-set gene in MAGE (x) and in the 358 Geuvadis "
  "individuals not in MAGE (y), at the same variant pair; colour = core or further main gene, shape = replication. Right: slope of "
  "each solid H1 example (lead variant in T2T-only sequence) in MAGE and Geuvadis.",
  "100 of 110 testable genes replicate (same direction, p < 0.05 with STAR and salmon), 86 after Bonferroni; core 53 of 54. "
  "8 of 15 solid H1 examples replicate beyond the best shared variant (5 distinct variants).",
  "The reference-driven changes are reproducible in independent individuals; H1 contributes a few solid examples."),
 ("Fig13 GWAS", "fig13-gwas.png", "Figure 13. GWAS lookup and colocalisation on each reference",
  "Left: colocalisation posterior PP4 (one shared causal variant) on the patched GRCh38 (x) and T2T (y) for gene x trait pairs "
  "with PP4 >= 0.5 on either; shaded = colocalises on one reference only. Right: share of genes with a GWAS Catalog hit in LD, "
  "and with a hit reached through only one reference's lead, for the main set and matched unchanged genes.",
  "No enrichment at GWAS loci (51% vs 46%), but the link depends on the reference far more often (26% vs 4%). 8 of 15 "
  "colocalising genes colocalise on one reference only; SuSiE-coloc confirms IFITM3-AS1, LOC124903630 and SERF1B, not PTGES3.",
  "The reference can change which gene a disease signal is attributed to."),
]
for name, png, title, how, what, take in figs:
    ws = wb.create_sheet(name)
    ws.column_dimensions["A"].width = 2; ws.column_dimensions["B"].width = 22; ws.column_dimensions["C"].width = 120
    ws.cell(1, 2, title).font = TITLE
    for i, (k, v) in enumerate([("How to read it", how), ("What it shows", what), ("Take-home", take)]):
        ws.cell(3 + i, 2, k).font = BOLD; c = ws.cell(3 + i, 3, v); c.font = BODY; c.alignment = WRAP
        ws.row_dimensions[3 + i].height = 30
    ws.cell(6, 3).fill = KEY
    ws.cell(7, 2, f"Source: results/07-presentation/figures/{png} (R, scripts/07-presentation/make-figures.R)").font = NOTE
    image(ws, "B9", FD / png, 1000)

# ------------------------------------------------------------------------------------------------------------
ws = wb.create_sheet("Key results")
ws.cell(1, 1, "Key results (MAGE discovery, 731 individuals; Geuvadis replication, 451; chr1-22)").font = TITLE
eg = pd.read_csv(DD / "egenes.tsv", sep="\t")
lab = {"chm13v2.0": "T2T-CHM13", "grch38": "GRCh38", "grch38-ffx": "Patched GRCh38"}
r = 3; ws.cell(r, 1, "eGenes (FDR 5%), matched variants").font = H2
rows = [[lab[x.reference], x.method, x.genes_tested, x.egenes, None] for x in eg.itertuples()]
r0 = r + 1
r = table(ws, r0, 1, ["Reference", "Counts", "Genes tested", "eGenes", "Share"], rows, widths=[26, 14, 14, 14, 14])
for i in range(len(rows)):
    c = ws.cell(r0 + 1 + i, 5, f"=D{r0+1+i}/C{r0+1+i}"); c.number_format = "0.0%"; c.font = BODY
    for j in (3, 4):
        ws.cell(r0 + 1 + i, j).number_format = "#,##0"
qcm = pd.read_csv(DD / "mapping-qc.tsv", sep="\t").filter(regex="pct").median()
ws.cell(r, 1, "Read mapping, median over 779 MAGE runs (STAR)").font = H2
rows = [[lab[k], round(qcm[f"unique_pct_{k}"], 2), round(qcm[f"multi_pct_{k}"], 2), round(qcm[f"unmapped_other_pct_{k}"], 2)] for k in lab]
r = table(ws, r + 1, 1, ["Reference", "Uniquely mapped %", "Multi-mapped %", "Unmapped 'other' %"], rows)
ws.cell(r - 1, 1, "Unmapped 'other' is higher on T2T (0.93% vs 0.09%); not yet explained (open item).").font = NOTE
r += 1
ws.cell(r, 1, "Genes changed between references (same gene, same variant; STAR and salmon agree; no artefact flag)").font = H2
rows = [["Patched GRCh38 vs GRCh38 (control)", 6, 7, 8, 0, 7],
        ["T2T vs GRCh38", 63, 121, 246, 59, 6],
        ["T2T vs patched GRCh38 (H0)", 57, 117, 255, 59, 0]]
r0 = r + 1
r = table(ws, r0, 1, ["Comparison", "Core (independent test)", "Main (paired, |delta| >= 0.1 SD)", "Extent (paired)",
                      "Main, in SDs", "Main, in FixItFelix regions"], rows)
ws.cell(r - 1, 1, "Sources: results/06-compare/eqtl-comparison/mage/gene-sets/changed-*.tsv (define-gene-sets.py); docs/t2t-eqtl.md 13B.4.").font = NOTE
r += 1
ws.cell(r, 1, "Paired test: calibration (300 random genes x 200 bootstrap resamples of individuals)").font = H2
cal = pd.read_csv(DD / "calibration.tsv", sep="\t")
rows = [[x.pair, x.method, round(x.median_boot_over_paired, 3), round(x.median_boot_over_indep, 3), round(x.lambda_indep, 3),
         round(x.lambda_paired, 2), round(x.median_rho, 3)] for x in cal.itertuples()]
r = table(ws, r + 1, 1, ["Pair", "Method", "Bootstrap SD / paired se", "Bootstrap SD / independent se", "lambda, independent",
                         "lambda, paired", "Median residual correlation"], rows)
ws.cell(r - 1, 1, "1 = the formula is right. The null pair (patched vs GRCh38, salmon) has lambda 1.03. Source: paired-matched/calibration.tsv.").font = NOTE
r += 1
ws.cell(r, 1, "How large a change? Gene count and replication by definition (T2T vs patched GRCh38)").font = H2
sf = pd.read_csv(DD / "size-floor.tsv", sep="\t")
rows = [[x.definition, x.genes, x.of_57, x.sd, x.tested, x.same_sign, x.nominal, x.bonferroni, None] for x in sf.itertuples()]
r0 = r + 1
r = table(ws, r0, 1, ["Definition", "Genes", "Of the 57", "In SDs", "Geuvadis tested", "Same sign", "p < 0.05", "Bonferroni",
                      "Share replicating"], rows)
for i in range(len(rows)):
    c = ws.cell(r0 + 1 + i, 9, f"=G{r0+1+i}/E{r0+1+i}"); c.number_format = "0%"; c.font = BODY
ws.cell(r - 1, 1, "Replication in the 358 Geuvadis individuals not in MAGE, paired test at MAGE's variant pairs, STAR and salmon. "
                  "Main set = |delta| >= 0.1 SD under both methods (decided 2026-10-29).").font = NOTE
r += 1
ws.cell(r, 1, "Replication in Geuvadis at the same variants (358 individuals not in MAGE)").font = H2
gs = pd.read_csv(DD / "gene-sets.tsv", sep="\t")
h0 = pd.read_csv(DD / "replication-h0.tsv", sep="\t").query("tested")
h1 = pd.read_csv(DD / "replication-h1.tsv", sep="\t")
h1s = h1[(h1.disc_genotype_call == "solid") & h1.disc_signal.isin(["T2T-only variant explains it", "both contribute"])]
rev = h0[h0.effect_class == "Reversed"]
rep = lambda x: [int(x.geuvadis_same_sign.notna().sum()), int(x.geuvadis_same_sign.eq(True).sum()),
                 int(x.geuvadis_nominal.eq(True).sum()), int(x.geuvadis_bonferroni.eq(True).sum())]
rows = [["Main set, 117 genes (paired test)"] + rep(gs[gs.tier.isin(["core", "main"])]),
        ["... core 57"] + rep(gs[gs.tier == "core"]),
        ["... further 60"] + rep(gs[gs.tier == "main"]),
        ["Extent, 255 genes"] + rep(gs),
        ["Reversals (core): slopes again of opposite sign, none near zero (STAR and salmon)", len(rev), int(reversed_again(rev).sum()), int(rev.nominal.sum()), int(rev.bonferroni.sum())],
        ["H1: solid T2T-only lead variants", len(h1s), int(h1s.same_sign.sum()), int(h1s.replicated.sum()), int(h1s.signal_replicated.sum())]]
r = table(ws, r + 1, 1, ["Set", "Tested", "Same direction", "Replicated (p < 0.05 both methods; H1: Bonferroni)",
                         "Bonferroni both methods; H1: also beyond shared variant"], rows)
ws.cell(r - 1, 1, "Sources: geuvadis/paired-matched-at-mage/ (paired-test.py --at mage) and replication-geuvadis-independent/ "
                  "(replicate-changes.py, replicate-h1.py --exclude-shared; the reversal row uses the independent test).").font = NOTE
r += 1
ws.cell(r, 1, "ctyper copy number as a covariate (genes covered by ctyper)").font = H2
cn = pd.read_csv(DD / "copy-number.tsv", sep="\t")
rows = [[name, len(x), int(x.still_changed.sum()), round(float((x.delta_cn / x.delta).median()), 2),
         ", ".join(x.loc[~x.still_changed, "symbol"]) or "—"]
        for name, x in [("Main set", cn[cn.tier.isin(["core", "main"])]), ("... core", cn[cn.tier == "core"]), ("Extent", cn)]]
r = table(ws, r + 1, 1, ["Set", "Covered", "Still changed, same sign", "Median delta ratio (with / without)", "Lost"], rows)
ws.cell(r - 1, 1, "Source: paired-matched-cn/ (paired-test.py --copy-number). NUTM2D stays nominally significant (STAR z -3.5).").font = NOTE
r += 1
ws.cell(r, 1, "Colocalisation with GWAS on each reference (gene x trait pairs with PP4 >= 0.5 on either)").font = H2
co = pd.read_csv(DD / "coloc.tsv", sep="\t")
su = pd.read_csv(ROOT / "results/06-compare/gwas/coloc-susie-mage-chm13v2.0-vs-grch38-ffx.tsv", sep="\t")
co = co.merge(su[["symbol", "trait", "susie_PP4_chm13v2.0", "susie_PP4_grch38-ffx"]], on=["symbol", "trait"], how="left")
co = co.sort_values(["call", "symbol"])
fmt = lambda v: None if pd.isna(v) else round(float(v), 3)
rows = [[x["symbol"], x["tier"], x["trait"], x["call"], fmt(x["PP4_chm13v2.0"]), fmt(x["PP4_grch38-ffx"]), fmt(x["PP3_chm13v2.0"]),
         fmt(x["PP3_grch38-ffx"]), x["gwas_min_p_chm13v2.0"], fmt(x["susie_PP4_chm13v2.0"]), fmt(x["susie_PP4_grch38-ffx"])]
        for _, x in co.iterrows()]
r0 = r + 1
r = table(ws, r0, 1, ["Gene", "Set", "Trait", "coloc call", "PP4 T2T", "PP4 patched", "PP3 T2T", "PP3 patched", "GWAS min p",
                      "SuSiE PP4 T2T", "SuSiE PP4 patched"], rows)
for i in range(len(rows)):
    ws.cell(r0 + 1 + i, 9).number_format = "0.0E+00"
ws.cell(r - 1, 1, "coloc.abf (one causal variant) and SuSiE-coloc (several; empty = no credible set for the eQTL or GWAS on that reference). "
                  "Sources: results/06-compare/gwas/coloc-*.tsv.").font = NOTE
r += 1
ws.cell(r, 1, "Positive control: q-value of the FixItFelix genes").font = H2
pcq = pd.read_csv(DD / "positive-control.tsv", sep="\t").pivot(index="symbol", columns="reference", values="qval")
rows = [[s, pcq.loc[s, "chm13v2.0"], pcq.loc[s, "grch38"], pcq.loc[s, "grch38-ffx"]] for s in pcq.sort_values("grch38").index]
r0 = r + 1
r = table(ws, r0, 1, ["Gene", "q, T2T", "q, GRCh38", "q, patched GRCh38"], rows)
for i in range(len(rows)):
    for j in (2, 3, 4):
        ws.cell(r0 + 1 + i, j).number_format = "0.0E+00"

# ------------------------------------------------------------------------------------------------------------
ws = wb.create_sheet("Gene sets")
ws.cell(1, 1, "Changed genes, T2T vs patched GRCh38: core (independent test) ⊂ main (paired, |delta| >= 0.1 SD) ⊂ extent (paired)").font = TITLE
gsx = pd.read_csv(DD / "gene-sets.tsv", sep="\t", dtype={"gene_ncbi": str})
gsx["tier_order"] = gsx.tier.map({"core": 0, "main": 1, "extent": 2})
gsx = gsx.sort_values(["tier_order", "z_paired"], key=lambda c: c.abs() if c.name == "z_paired" else c, ascending=[True, False])
cls = {"sd": "Segmental duplication", "other": "Other", "t2t_new_sequence": "T2T-only sequence", "ffx_patched": "FixItFelix region"}
yn = lambda v: None if pd.isna(v) else ("yes" if v in (True, "True") else "no")
rows = [[x.symbol, x.gene_ncbi, x.biotype, cls.get(x.region_class, x.region_class), x.tier, round(x.slope_a, 3), round(x.slope_b, 3),
         round(x.delta, 3), round(x.z_paired, 1), round(x.delta_salmon, 3), round(x.z_paired_salmon, 1),
         yn(x.geuvadis_nominal), yn(x.geuvadis_bonferroni)] for x in gsx.itertuples()]
r = table(ws, 3, 1, ["Gene", "NCBI GeneID", "Biotype", "Region class", "Set", "Slope T2T", "Slope patched GRCh38", "Delta (STAR)",
                     "Paired z (STAR)", "Delta (salmon)", "Paired z (salmon)", "Geuvadis p < 0.05", "Geuvadis Bonferroni"],
          rows, widths=[16, 12, 20, 20, 9, 10, 12, 11, 11, 11, 11, 12, 12])
ws.auto_filter.ref = f"A3:M{3+len(rows)}"
ws.cell(r, 1, "Set 'main' = further genes of the main set (the main set is core + main, 117 genes); 'extent' = changed but |delta| < 0.1 SD "
              "under at least one method. Geuvadis: 358 individuals not in MAGE, MAGE's variant pair, paired test, STAR and salmon "
              "(empty = not testable). Source: results/06-compare/eqtl-comparison/mage/gene-sets/changed-chm13v2.0-vs-grch38-ffx.tsv.").font = NOTE

ws = wb.create_sheet("57 genes")
ws.cell(1, 1, "The 57 core genes: changed with the independent-samples test (STAR and salmon agree, no artefact flag)").font = TITLE
rg = pd.read_csv(DD / "robust-genes.tsv", sep="\t", dtype={"gene_ncbi": str}).sort_values("p_diff")
cls = {"sd": "Segmental duplication", "other": "Other", "t2t_new_sequence": "T2T-only sequence", "ffx_patched": "FixItFelix region"}
rows = [[x.symbol, x.gene_ncbi, x.biotype, cls.get(x.region_class, x.region_class),
         {"stronger on chm13v2.0": "Stronger on T2T", "stronger on grch38-ffx": "Stronger on patched GRCh38"}.get(x.pattern, "Opposite sign"), round(x.slope_a, 3), round(x.slope_b, 3),
         round(x.z_diff, 1), round(x.z_diff_salmon, 1), x.p_diff, round(x.geno_r, 3),
         round(x.cpm_chm13v2_0, 2) if hasattr(x, "cpm_chm13v2_0") else None, None,
         x.partner, cls.get(x.partner_class, x.partner_class), round(x.partner_r, 2) if pd.notna(x.partner_r) else None]
        for x in rg.rename(columns={"cpm_chm13v2.0": "cpm_chm13v2_0", "cpm_grch38-ffx": "cpm_grch38_ffx"}).itertuples()]
cpm_b = rg["cpm_grch38-ffx"].round(2).tolist()
for row, b in zip(rows, cpm_b):
    row[12] = b
hdr = ["Gene", "NCBI GeneID", "Biotype", "Region class", "Effect pattern", "Slope T2T", "Slope patched GRCh38", "z (STAR)",
       "z (salmon)", "p (STAR)", "Genotype r between panels", "CPM T2T", "CPM patched GRCh38", "Read partner",
       "Partner class", "Partner r"]
r = table(ws, 3, 1, hdr, rows, widths=[16, 12, 20, 20, 15, 10, 12, 9, 9, 10, 12, 10, 12, 18, 20, 9])
for i in range(len(rows)):
    ws.cell(4 + i, 10).number_format = "0.0E+00"
ws.auto_filter.ref = f"A3:P{3+len(rows)}"
ws.cell(r, 1, "Slopes are the effect of the same variant (alternative allele, T2T coding) on normalised expression. "
              "Read partner = the gene whose per-person change mirrors this gene's (most negative correlation); "
              "r < -0.5 counts as a clear partner. Genotype r of -1 means the two panels code the alleles the other way round (same variant). Sources: robust-matched/changed-*.tsv and read-destination-*.tsv.").font = NOTE

# ------------------------------------------------------------------------------------------------------------
ws = wb.create_sheet("Glossary")
ws.cell(1, 1, "Glossary").font = TITLE
gl = [("eQTL", "A DNA variant whose alleles go with different expression levels of a nearby gene (cis: within 1 Mb)."),
      ("eGene", "A gene with at least one significant eQTL (here FDR 5% after permutations)."),
      ("Lead variant", "The variant with the strongest association for a gene."),
      ("Slope / effect size", "Change in normalised expression per copy of the alternative allele."),
      ("z (effect difference)", "(slope on reference A - slope on B) / standard error of the difference, for the same gene and variant."),
      ("Paired test", "The effect difference tested with the covariance of the two slopes, because both come from the same individuals; "
       "the independent-samples test ignores it and is far too conservative here."),
      ("Gene sets: core / main / extent", "Core: 57 genes changed with the old independent test. Main: 117 genes changed with the paired "
       "test by at least 0.1 SD under both methods (includes core). Extent: 255 genes changed with the paired test, any size."),
      ("Bootstrap", "Resampling individuals with replacement to measure the spread of an estimate directly."),
      ("Genomic inflation (lambda)", "Median test statistic over all genes relative to its expected value under no change; 1 = as expected."),
      ("Copy number (ctyper)", "How many copies of a gene (or gene family) a person carries, genotyped from short reads with pangenome alleles."),
      ("GWAS", "Genome-wide association study: variants associated with a disease or trait."),
      ("Colocalisation (coloc)", "Test whether an eQTL and a GWAS signal share one causal variant. PP4 = posterior probability of a "
       "shared variant; PP3 = two different variants."),
      ("SuSiE-coloc", "Colocalisation allowing several causal variants per region, using LD between variants."),
      ("LD (r2)", "Linkage disequilibrium: how strongly the genotypes of two variants go together in a population."),
      ("Allele-specific expression", "Within one heterozygous person, more RNA from one haplotype than the other; tests an eQTL's direction."),
      ("LCL", "Lymphoblastoid cell line: immortalised B cells, the cell type of MAGE, Geuvadis and GM12878."),
      ("T2T-CHM13", "The first complete (telomere-to-telomere) human genome assembly (Nurk 2022)."),
      ("GRCh38", "The standard human reference since 2013; has gaps, false duplications and collapsed duplications."),
      ("FixItFelix / patched GRCh38", "GRCh38 with ~20 known errors fixed (false copies masked, missing paralogs added as decoys; Behera 2023)."),
      ("Segmental duplication (SD)", "A long stretch (> 1 kb) duplicated elsewhere in the genome with high identity (here >= 98%)."),
      ("Paralog", "A gene copy created by duplication; reads from one copy can map to another."),
      ("Processed pseudogene", "A copy of a gene's mRNA re-inserted into the genome; it attracts reads from the parent gene."),
      ("Multi-mapped read", "A read that aligns equally well to several places."),
      ("STAR / salmon", "Two ways to count reads per gene: STAR counts unique alignments per gene; salmon distributes reads across transcripts."),
      ("Matched variant", "A variant called independently on both references and paired by position and alleles."),
      ("Robust change", "Changed with STAR and with salmon, same direction, no artefact flag (annotation overlap, MHC, PAR)."),
      ("Read destination / partner", "The gene that gains the reads a changed gene loses (or vice versa), found by correlating per-person changes."),
      ("FDR / q-value", "False discovery rate; q <= 0.05 means about 5% of the calls are expected to be false."),
      ("CPM", "Counts per million reads.")]
table(ws, 3, 1, ["Term", "Meaning"], gl, widths=[30, 120])

# ------------------------------------------------------------------------------------------------------------
ws = wb.create_sheet("Status & next")
ws.cell(1, 1, "Status, caveats and next steps").font = TITLE
st = [("Done", "MAGE expression on three references (779 runs, 731 individuals)", "2026-10-24"),
      ("Done", "cis-eQTL mapping on MAGE, chr1-22, STAR and salmon, matched variants", "2026-10-25"),
      ("Done", "Same-variant comparison, robust set (57 genes), read destination (33 with a clear partner)", "2026-10-25"),
      ("Done", "eQTL mapping with all variants; H1 lead-variant quality check (15 solid examples)", "2026-10-26"),
      ("Done", "Geuvadis expression (652 runs, 451 individuals) and eQTL mapping on three references, STAR and salmon", "2026-10-27"),
      ("Done", "Replication at the same variants: 50/54 changed genes (44 Bonferroni), 10/12 reversals, 8 H1 genes at 5 variants", "2026-10-28"),
      ("Done", "Replication without the 93 individuals in both cohorts (358): unchanged", "2026-10-28"),
      ("Done", "Paired test of effect changes, bootstrap calibration; gene sets core 57 / main 117 / extent 255", "2026-10-28/29"),
      ("Done", "ctyper copy number as covariate: 55 of 62 covered main-set genes unchanged", "2026-10-28"),
      ("Done", "GWAS Catalog lookup with matched background; coloc with six GWAS on both references; SuSiE-coloc", "2026-10-29"),
      ("Running", "Allele-specific expression: RNA realignment of 347 individuals (target genes) and ONT SNP calls in 260", "2026-10-29"),
      ("Next", "Allele-specific counts, reference-bias correction and direction test per gene and reference", ""),
      ("Next", "Read destination for the 60 further main genes, and in Geuvadis", ""),
      ("Next", "chrX (pseudo-autosomal boundaries per reference)", ""),
      ("Next", "Genes without a single read partner (24 core): look at multi-mapped reads in the BAMs", ""),
      ("Next", "Link changed eQTLs to GM12878 regulatory elements (t2t-atlas, H5)", ""),
      ("Caveat", "salmon alone calls ~550 genes changed, mostly not seen with STAR: transcript annotations differ between the two assemblies; hence the both-methods rule", ""),
      ("Caveat", "An allele-orientation error for some indels (fixed 2026-10-25) had added 12 false genes; orientation now comes from the genotypes", ""),
      ("Caveat", "'Stronger on' (larger |slope|) is not meaningful when the slopes have opposite signs (31 of 57 genes); tables now give the effect pattern", ""),
      ("Caveat", "Single-variant coloc is not robust at PTGES3 and SCAF11 (not confirmed by SuSiE-coloc); 4 loci untestable with European reference LD", ""),
      ("Caveat", "Copy number explains part of the change for NUTM2D, PDPR2P, SMN1, SERF1B, NPIPB6, GRAPL, GSTM2", ""),
      ("Caveat", "Unmapped 'other' reads are higher on T2T (0.93% vs 0.09%); not yet explained", ""),
      ("Caveat", "LOC124906110 correlates with many unrelated genes; treat as a possible technical artefact", "")]
table(ws, 3, 1, ["Status", "Item", "Date"], st, widths=[12, 120, 20])

for s in wb.worksheets:
    s.sheet_view.showGridLines = False
# no cached formula values are written (openpyxl), so Excel recalculates everything on open
from openpyxl.workbook.properties import CalcProperties
wb.calculation = CalcProperties(fullCalcOnLoad=True)
OUT.parent.mkdir(parents=True, exist_ok=True)
wb.save(OUT)
print(OUT)
