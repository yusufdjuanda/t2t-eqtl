"""Presentation: one Excel workbook that explains the project, its data and its first results (oral-talk format).

    pixi run python scripts/07-presentation/build-workbook.py

Inputs: results/07-presentation/data/*.tsv (prepare-figure-data.py), results/07-presentation/figures/*.png
(pipeline.dot via Graphviz; make-figures.R). Output: results/07-presentation/t2t-eqtl-project-overview.xlsx.
Text in the sheets follows docs/t2t-eqtl.md (sections 1-4, 12, 13B.3); update both together.
"""
from datetime import date
from pathlib import Path

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
          "4. The patch changes 8 genes, all in its regions; T2T changes 57 more genes outside them (robust to the "
          "quantification method), mostly by moving reads between paralogs, pseudogenes and genes that exist only in T2T."]:
    r = text(ws, r, 2, s, BODY, 3)
r += 1
r = text(ws, r, 2, "Key numbers", H2)
key = [("Individuals (discovery, MAGE)", 731), ("Individuals (replication, Geuvadis)", 451),
       ("References compared", "T2T-CHM13v2.0, GRCh38 no-alt, FixItFelix patched GRCh38"),
       ("eGenes per reference (FDR 5%, STAR)", "15,013 / 15,259 / 15,224"),
       ("Genes changed by the patch (control)", 8), ("Genes changed by T2T beyond the patch (robust)", 57),
       ("... of which with a clear read destination", 33)]
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
          ("Fig2 ... Fig10", "One sheet per figure (R, ggplot2) with a caption: what it shows, how to read it, the take-home"),
          ("Key results", "The numbers behind the figures"),
          ("57 genes", "The robust gene set with effect sizes and where each gene's reads go"),
          ("Glossary", "Terms used in the talk, in plain language"),
          ("Status & next", "What is done, running and next; known caveats")]
r = table(ws, r, 2, ["Sheet", "What is in it"], sheets)

# ------------------------------------------------------------------------------------------------------------
ws = wb.create_sheet("Talk outline")
ws.cell(1, 1, "Oral presentation plan (about 15 minutes + questions)").font = TITLE
slides = [
 (1, "Title", "A complete reference genome changes the eQTL map beyond known GRCh38 fixes", "—", 0.5,
  "Name, affiliation. One line: same raw RNA-seq reads, three references, what changes."),
 (2, "Why the reference matters", "Every step of an eQTL study (genotypes, read placement, gene models) depends on the reference genome",
  "—", 1.0, "GRCh38 has 151 Mb of unknown sequence, false duplications (chr21) and collapsed duplications. T2T-CHM13 adds or "
  "corrects 238 Mb (Nurk 2022). Reads from missing sequence do not vanish: they land on the most similar place that exists."),
 (3, "What is already known", "FixItFelix patched ~20 GRCh38 errors and found artefactual eQTLs (PWP2, GATD3A) in Geuvadis",
  "—", 1.0, "Behera et al. 2023 Genome Biology. Limits: patched regions only (a few Mb), GRCh38 genotypes and gene models. "
  "No eQTL study has been run natively on T2T (searched Oct 2026)."),
 (4, "Question", "H0: most reference-driven eQTL changes lie OUTSIDE the regions FixItFelix patched", "—", 0.5,
  "The patched GRCh38 is our third reference: it separates what known patches fix from what a complete genome adds."),
 (5, "Data", "731 MAGE LCLs from 26 populations (discovery) and 451 Geuvadis LCLs (replication), all with 1000 Genomes genotypes",
  "Fig2", 1.0, "MAGE: one lab, 2x150 bp, ~28 M read pairs per run. Geuvadis: several labs, 2x76 bp. 93 people are in both. "
  "Genotypes: the phased 1KGP panel called natively on T2T, and the NYGC panel on GRCh38."),
 (6, "Design", "Same raw reads processed three times; only the reference changes", "Fig1", 1.5,
  "SRA download with checksum and read-count checks; STAR + salmon per reference; tensorQTL with 10,000 permutations, "
  "5 genotype PCs, sex and 60 expression PCs. Comparison on variants present on both references (matched)."),
 (7, "Mapping", "T2T halves multi-mapping reads (4.2% -> 2.1%) and places slightly more reads uniquely", "Fig3", 1.0,
  "Grey lines join the same run on each reference. Fewer ambiguous reads is the first sign that read placement changes."),
 (8, "Fair comparison of genotypes", "12.6 M variants are matched between the T2T and GRCh38 panels; the main test uses only these",
  "Fig4", 1.0, "LiftOver is used only to pair variants called independently on each reference. Matched pairs must agree "
  "in allele frequency and in person-level genotypes, so any change comes from expression, not genotypes."),
 (9, "Overall", "About 15,000 eGenes on every reference: the overall map is stable", "Fig5", 0.5,
  "Agrees with the MAGE paper (~15,000 on GRCh38). The interesting part is which genes differ."),
 (10, "Positive control", "The known chr21 artefacts reappear: PWP2, GATD3, CBS, ICOSLG, U2AF1 are inflated only on unpatched GRCh38",
  "Fig6", 1.0, "T2T agrees with the patch. The patch changes 9 genes (8 with both quantifiers), all in its regions. "
  "This independently reproduces FixItFelix in a larger, more diverse cohort."),
 (11, "Genome-wide comparison", "Same gene, same variant: 57 genes change between T2T and the patched GRCh38, none in FixItFelix regions",
  "Fig7", 1.5, "Robust = changed with STAR AND salmon, same direction, no annotation/MHC/PAR flag. Effects are compared "
  "on the same variant; the allele orientation comes from the genotypes themselves."),
 (12, "The 57 genes", "Mostly segmental duplications and gene families (CCZ1, ASAH2, PMS2, NBPF, NPIP, TBC1D3, STAG3L) plus "
  "highly expressed genes with pseudogenes (RPL9, RPSA, UBC)", "Fig8", 1.0,
  "About 8 times as many genes as the patch fixes. Half have opposite-sign effects at the same variant."),
 (13, "Mechanism", "One gene loses the reads its partner gains: paralog swaps, processed pseudogenes and genes that exist only in T2T",
  "Fig9", 1.5, "Per person, CCZ1 and CCZ1B change in mirror image (r = -0.96). PDPR on GRCh38 holds part of the reads of a "
  "gene in T2T-only sequence (r = -0.82, about 8% of its change): missing sequence creates an eQTL artefact outside every patch."),
 (14, "What it means for an eQTL", "The same people and the same variant give a different eQTL depending on the reference",
  "Fig10", 1.0, "PWP2: the known artefact. CCZ1: the effect reverses. SMN2: much stronger on T2T. These genes are tested "
  "in GTEx and used in colocalisation; the reference choice changes the answer."),
 (15, "Limits and next steps", "Replication in Geuvadis (running); T2T-only variants; chrX; paralog copy number (ctyper)",
  "—", 1.0, "Salmon alone over-calls changes (transcript annotation differs between assemblies), so both quantifiers are "
  "required. 24 of 57 genes have no single read partner yet."),
 (16, "Conclusions", "A complete reference changes eQTLs well beyond known patches, through read placement in duplicated "
  "and newly resolved sequence", "—", 0.5, "Take-home: use T2T (or at least check SD genes) before trusting an eQTL in "
  "a gene family. Code and results will be released."),
]
r = table(ws, 3, 1, ["Slide", "Title", "Key message (one sentence on the slide)", "Figure", "Minutes", "Speaker notes (plain language)"],
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
  "scripts/03-expression/fetch-sra-batch.sbatch (prefetch)", "vdb-validate (NCBI checksums) + read pairs equal to ENA", "Expression on all three references"),
 ("RNA-seq", "Geuvadis RNA-seq (replication cohort)", "LCL RNA-seq of 462 individuals, 5 populations, 7 labs; 652 runs used (451 individuals with T2T genotypes), 2x76 bp, unstranded",
  "Lappalainen et al. 2013, Nature", "ENA PRJEB3366 (mirrored at NCBI SRA); sample sheet E-GEUV-1 (BioStudies)", "2026-10-17", "~2.5 TB reads",
  "scripts/03-expression/fetch-sra-batch.sbatch", "vdb-validate + read pairs equal to ENA", "Replication (running since 2026-10-25)"),
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
  "t2t-atlas prepare-references.sbatch; scripts/02-rnaseq-pilot/prepare-star-references.sbatch", "CHECKSUMS.md5", "Reference 1"),
 ("Reference", "GRCh38 no-alt analysis set", "Primary assembly + unplaced contigs, rCRS chrM, chrEBV, chrY PARs masked (the set ENCODE uses)",
  "Genome Reference Consortium / NCBI", "NCBI GCA_000001405.15 seqs_for_alignment_pipelines", "2026-10-04", "3.1 Gb",
  "t2t-atlas prepare-noalt-reference.sbatch", "NCBI md5checksums.txt", "Reference 2"),
 ("Reference", "FixItFelix patched GRCh38", "GRCh38 no-alt with false duplications masked and decoys for collapsed genes",
  "Behera et al. 2023, Genome Biology", "Zenodo 7532965 (GRCh38_masked_v2_decoy_gene.fasta.gz)", "2026-10-17", "873 MB",
  "scripts/03-expression/prepare-ffx-reference.sbatch", "md5 35d00c3c48d1849adb8e02086c3a1c70", "Reference 3 (separates patch effects from completeness)"),
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
  "curl", "md5 equal to Zenodo", "Planned covariate (copy number)"),
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
ws.cell(note_r, 1, "Mapping rates (STAR) are filled for MAGE; Geuvadis runs are being processed. Two Geuvadis runs "
                   "(NA07000, HG00237) have no population in E-GEUV-1. Source: config/samplesheets/*.csv, "
                   "results/03-expression/mage/matrices/runs.tsv.").font = NOTE

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
 ("2 / 2b", "Pilots", "5 Geuvadis + 5 MAGE runs on both references; go/no-go 1 passed", "STAR, salmon", "scripts/02-rnaseq-pilot/", "~2 h"),
 ("3", "Expression (Nextflow)", "Per run: SRA download and checks, STAR GeneCounts and salmon on three references, read counts per contig; "
  "merged into gene x individual matrices", "sra-tools, STAR 2.7.11b, salmon 1.11.4, Nextflow", "scripts/03-expression/ (nextflow/, submit-cohort.sh)",
  "MAGE: 6,157 CPU h (779 runs x 3 references)"),
 ("4", "Genotypes on GRCh38", "Download and subset the NYGC phased panel", "curl, bcftools", "scripts/01-genotypes/grch38/", "hours"),
 ("6a", "Variant matching", "Pair T2T and GRCh38 variants (liftOver only to pair, same alleles, MAF within 0.02); 12.6 M matched", "UCSC liftOver, pandas",
  "scripts/06-compare/match-variants.*", "~1 h"),
 ("5", "cis-eQTL mapping", "TMM + inverse-normal expression; 5 genotype PCs, sex, 60 expression PCs; map_cis (10,000 permutations) and map_nominal per chromosome",
  "tensorQTL 1.0.10, plink2", "scripts/05-eqtl/", "~44 CPU h per reference and method"),
 ("6", "Comparison", "eGene presence, same-variant effect sizes (z), STAR AND salmon rule, region classes and artefact flags, read destination",
  "pandas, scipy", "scripts/06-compare/", "minutes"),
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
  "Each point is a gene: the effect (slope) of its lead variant on patched GRCh38 (x) and on T2T (y). On the diagonal = no change.",
  "Almost all of ~15,000 genes lie on the diagonal. 57 change robustly (coloured): with both quantifiers, same direction, no "
  "artefact flag. None lies in a FixItFelix region; half are in segmental duplications.",
  "T2T changes eQTLs beyond the known patches (H0)."),
 ("Fig8 57 genes", "fig08-robust-genes.png", "Figure 8. The 57 robustly changed genes",
  "Left: effect difference z per gene (colour = region class). Right top: genes by which reference shows the stronger effect. "
  "Right bottom: how strongly each gene's read loss is mirrored by one partner gene.",
  "Changes go both ways (30 stronger on T2T, 27 on patched GRCh38). 33/57 have a clear partner gene that takes their reads.",
  "The changes are a read-placement phenomenon in gene families."),
 ("Fig9 Read swaps", "fig09-paralog-swaps.png", "Figure 9. Where the reads go",
  "Each dot is one person: change in the gene's expression (T2T minus patched GRCh38) vs change in its partner.",
  "CCZ1/CCZ1B and ASAH2/ASAH2B swap reads almost exactly; RPSA loses reads to its processed pseudogene RPSA2; "
  "PDPR on GRCh38 holds part of the reads that T2T assigns to a gene in sequence GRCh38 lacks (r = -0.82, about 8% of its change).",
  "Missing or collapsed sequence moves reads between copies, and the move differs between people."),
 ("Fig10 Examples", "fig10-eqtl-examples.png", "Figure 10. Three eQTLs on two references",
  "Normalised expression by genotype at the lead variant, for the same people, on two references.",
  "PWP2: strong genotype effect only on unpatched GRCh38 (the known artefact). CCZ1: the effect reverses direction. "
  "SMN2: a much stronger effect on T2T.",
  "The reference can create, remove or reverse an eQTL."),
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
ws.cell(1, 1, "Key results (MAGE, 731 individuals, chr1-22)").font = TITLE
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
ws.cell(r, 1, "Genes changed between references (same gene, same variant; Bonferroni 5%)").font = H2
rows = [["Patched GRCh38 vs GRCh38 (control)", 9, 9, 8, 6, 0, 6],
        ["T2T vs GRCh38", 145, 563, 79, 63, 34, 5],
        ["T2T vs patched GRCh38 (H0)", 138, 549, 71, 57, 33, 0]]
r0 = r + 1
r = table(ws, r0, 1, ["Comparison", "STAR", "salmon", "STAR and salmon", "... no artefact flag", "... in SDs", "... in FixItFelix regions"], rows)
for j, w in zip(range(6, 8), (12, 20)):
    ws.column_dimensions[get_column_letter(j)].width = w
ws.cell(r - 1, 1, "Source: results/06-compare/eqtl-comparison/mage/robust-matched/summary.tsv (after the allele-sign fix of 2026-10-25).").font = NOTE
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
ws = wb.create_sheet("57 genes")
ws.cell(1, 1, "The robust set: genes T2T changes relative to the patched GRCh38 (STAR and salmon agree, no artefact flag)").font = TITLE
rg = pd.read_csv(DD / "robust-genes.tsv", sep="\t", dtype={"gene_ncbi": str}).sort_values("p_diff")
cls = {"sd": "Segmental duplication", "other": "Other", "t2t_new_sequence": "T2T-only sequence", "ffx_patched": "FixItFelix region"}
rows = [[x.symbol, x.gene_ncbi, x.biotype, cls.get(x.region_class, x.region_class),
         "T2T" if x.stronger_on == "chm13v2.0" else "Patched GRCh38", round(x.slope_a, 3), round(x.slope_b, 3),
         round(x.z_diff, 1), round(x.z_diff_salmon, 1), x.p_diff, round(x.geno_r, 3),
         round(x.cpm_chm13v2_0, 2) if hasattr(x, "cpm_chm13v2_0") else None, None,
         x.partner, cls.get(x.partner_class, x.partner_class), round(x.partner_r, 2) if pd.notna(x.partner_r) else None]
        for x in rg.rename(columns={"cpm_chm13v2.0": "cpm_chm13v2_0", "cpm_grch38-ffx": "cpm_grch38_ffx"}).itertuples()]
cpm_b = rg["cpm_grch38-ffx"].round(2).tolist()
for row, b in zip(rows, cpm_b):
    row[12] = b
hdr = ["Gene", "NCBI GeneID", "Biotype", "Region class", "Stronger on", "Slope T2T", "Slope patched GRCh38", "z (STAR)",
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
      ("z (effect difference)", "(slope on reference A - slope on B) / combined standard error, for the same gene and variant."),
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
      ("Running", "Geuvadis expression (652 runs, 11 batches) for replication", "expected 2026-10-27"),
      ("Running", "eQTL mapping with all variants, including T2T-only variants (H1)", "2026-10-25"),
      ("Next", "chrX (pseudo-autosomal boundaries per reference)", ""),
      ("Next", "Paralog copy number from ctyper as a covariate", ""),
      ("Next", "Genes without a single read partner (24): look at multi-mapped reads in the BAMs", ""),
      ("Next", "Link changed eQTLs to GM12878 regulatory elements (t2t-atlas) and GWAS colocalisation (H5, H6)", ""),
      ("Caveat", "salmon alone calls ~550 genes changed, mostly not seen with STAR: transcript annotations differ between the two assemblies; hence the both-methods rule", ""),
      ("Caveat", "An allele-orientation error for some indels (fixed 2026-10-25) had added 12 false genes; orientation now comes from the genotypes", ""),
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
