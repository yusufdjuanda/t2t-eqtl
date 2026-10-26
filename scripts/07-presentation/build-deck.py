"""Presentation: the oral-talk deck (PowerPoint) built from the figures and the workbook's content.

    pixi run python scripts/07-presentation/build-deck.py

Inputs: results/07-presentation/figures/*.png (pipeline.dot, make-figures.R) and results/07-presentation/data/*.tsv.
Output: results/07-presentation/t2t-eqtl-talk.pptx (16:9, 18 slides, speaker notes on every slide).
Text follows build-workbook.py (Talk outline) and docs/t2t-eqtl.md; update them together.
"""
from pathlib import Path

import pandas as pd
from PIL import Image as PILImage
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.oxml.ns import qn
from pptx.util import Emu, Inches, Pt

ROOT = Path(__file__).resolve().parents[2]
FD = ROOT / "results/07-presentation/figures"
DD = ROOT / "results/07-presentation/data"
OUT = ROOT / "results/07-presentation/t2t-eqtl-talk.pptx"

HEAD, BODY = "Cambria", "Calibri"
NAVY, INK, MUTED, TINT, WHITE = "14213D", "1F2933", "5B6770", "EEF2F7", "FFFFFF"
T2T, G38, FFX, SD, NEW = "0072B2", "D55E00", "E69F00", "CC79A7", "009E73"
ICE = "C9D6EA"
rgb = lambda h: RGBColor.from_string(h)

prs = Presentation()
prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
BLANK = prs.slide_layouts[6]
W, H = 13.333, 7.5
state = {"n": 0}


def slide(dark=False):
    s = prs.slides.add_slide(BLANK)
    s.background.fill.solid(); s.background.fill.fore_color.rgb = rgb(NAVY if dark else WHITE)
    state["n"] += 1
    return s


def box(s, x, y, w, h, name=None):
    tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame; tf.word_wrap = True; tf.auto_size = None
    tf.margin_left = tf.margin_right = Inches(0.05); tf.margin_top = tf.margin_bottom = Inches(0.03)
    if name:
        tb.name = name
    return tb, tf


def run(p, text, size, color=INK, bold=False, italic=False, font=BODY):
    r = p.add_run(); r.text = text
    f = r.font; f.size = Pt(size); f.bold = bold; f.italic = italic; f.name = font; f.color.rgb = rgb(color)
    return r


def para(tf, first=False, align=PP_ALIGN.LEFT, after=0):
    p = tf.paragraphs[0] if first else tf.add_paragraph()
    p.alignment = align; p.space_after = Pt(after)
    return p


def text(s, x, y, w, h, parts, size=16, color=INK, bold=False, align=PP_ALIGN.LEFT, font=BODY, anchor=MSO_ANCHOR.TOP,
         italic=False, name=None):
    """parts: a string, or a list of (text, {overrides}) runs in one paragraph."""
    tb, tf = box(s, x, y, w, h, name); tf.vertical_anchor = anchor
    p = para(tf, True, align)
    for t, o in ([(parts, {})] if isinstance(parts, str) else parts):
        run(p, t, o.get("size", size), o.get("color", color), o.get("bold", bold), o.get("italic", italic), o.get("font", font))
    return tb


def bullets(s, x, y, w, h, items, size=16, color=INK, after=8, name=None):
    """items: strings or (bold lead, rest) tuples; real bullets via buChar."""
    tb, tf = box(s, x, y, w, h, name)
    for i, it in enumerate(items):
        p = para(tf, i == 0, after=after)
        pPr = p._p.get_or_add_pPr(); pPr.set("marL", str(Inches(0.28))); pPr.set("indent", str(-Inches(0.24)))
        bu = pPr.makeelement(qn("a:buChar"), {"char": "•"}); pPr.append(bu)
        if isinstance(it, tuple):
            run(p, it[0], size, color, bold=True); run(p, it[1], size, color)
        else:
            run(p, it, size, color)
    return tb


def title(s, t, kicker=None, dark=False):
    text(s, 0.6, 0.35, 12.1, 0.8, t, size=32, bold=True, font=HEAD, color=WHITE if dark else NAVY, name="Title")
    if kicker:
        text(s, 0.6, 1.1, 12.1, 0.5, kicker, size=17, italic=True, color=ICE if dark else MUTED, name="Kicker")


def footer(s):
    text(s, 0.6, 7.05, 8, 0.3, "T2T-native eQTL study  ·  MAGE, 731 lymphoblastoid cell lines", size=10, color=MUTED, name="Footer")
    text(s, 11.9, 7.05, 0.85, 0.3, str(state["n"]), size=10, color=MUTED, align=PP_ALIGN.RIGHT, name="Slide number")


def card(s, x, y, w, h, fill=TINT, name=None):
    sh = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(x), Inches(y), Inches(w), Inches(h))
    sh.adjustments[0] = 0.06
    sh.fill.solid(); sh.fill.fore_color.rgb = rgb(fill); sh.line.fill.background(); sh.shadow.inherit = False
    if name:
        sh.name = name
    return sh


def dot(s, x, y, d, fill, label=None, size=14, name=None):
    sh = s.shapes.add_shape(MSO_SHAPE.OVAL, Inches(x), Inches(y), Inches(d), Inches(d))
    sh.fill.solid(); sh.fill.fore_color.rgb = rgb(fill); sh.line.fill.background(); sh.shadow.inherit = False
    if label:
        tf = sh.text_frame; tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
        tf.vertical_anchor = MSO_ANCHOR.MIDDLE
        p = tf.paragraphs[0]; p.alignment = PP_ALIGN.CENTER; run(p, label, size, WHITE, bold=True)
    if name:
        sh.name = name
    return sh


def stat(s, x, y, w, big, label, color=T2T, big_size=44, label_size=14, dark=False):
    bh = max(0.85, big_size / 72 * 1.3)
    text(s, x, y, w, bh, big, size=big_size, bold=True, color=color, font=HEAD)
    text(s, x, y + bh, w, 0.75, label, size=label_size, color=ICE if dark else MUTED)


def picture(s, png, x, y, w, h, align="center", name=None):
    """Fit png inside the box keeping its aspect; returns the placed (x, y, w, h) in inches."""
    iw, ih = PILImage.open(png).size
    scale = min(w / iw, h / ih); pw, ph = iw * scale, ih * scale
    px = x + (w - pw) / 2 if align == "center" else x
    py = y + (h - ph) / 2
    pic = s.shapes.add_picture(str(png), Inches(px), Inches(py), Inches(pw), Inches(ph))
    pic.name = name or png.stem
    pic._element.nvPicPr.cNvPr.set("descr", name or png.stem)
    return px, py, pw, ph


def table(s, x, y, w, rows, col_w, size=12, head_fill=NAVY, row_h=0.42, name=None):
    gt = s.shapes.add_table(len(rows), len(rows[0]), Inches(x), Inches(y), Inches(w), Inches(row_h * len(rows)))
    if name:
        gt.name = name
    t = gt.table
    t.first_row = True; t.horz_banding = False
    for j, cw in enumerate(col_w):
        t.columns[j].width = Inches(cw)
    for i, r in enumerate(rows):
        t.rows[i].height = Inches(row_h)
        for j, v in enumerate(r):
            c = t.cell(i, j); c.text = ""
            c.margin_left = c.margin_right = Inches(0.08); c.margin_top = c.margin_bottom = Inches(0.03)
            c.vertical_anchor = MSO_ANCHOR.MIDDLE
            c.fill.solid(); c.fill.fore_color.rgb = rgb(head_fill if i == 0 else (TINT if i % 2 == 0 else WHITE))
            p = c.text_frame.paragraphs[0]
            run(p, str(v), size, WHITE if i == 0 else INK, bold=(i == 0 or j == 0))
    return t


def notes(s, t):
    s.notes_slide.notes_text_frame.text = t


fig = lambda n: next(FD.glob(f"fig{n:02d}-*.png"))

# 1 ---------------------------------------------------------------------------------------------------------
s = slide(dark=True)
text(s, 0.8, 1.4, 11.5, 2.3, "A complete reference genome changes the eQTL map beyond known GRCh38 fixes",
     size=40, bold=True, font=HEAD, color=WHITE, name="Title")
text(s, 0.8, 3.85, 11.5, 0.6, "T2T-native cis-eQTL mapping in 731 MAGE lymphoblastoid cell lines", size=20, color=ICE)
for i, (c, lab) in enumerate([(T2T, "T2T-CHM13"), (G38, "GRCh38"), (FFX, "Patched GRCh38")]):
    dot(s, 0.8 + i * 2.6, 4.85, 0.32, c)
    text(s, 1.2 + i * 2.6, 4.82, 2.2, 0.4, lab, size=15, color=WHITE)
text(s, 0.8, 5.35, 11, 0.4, "Same raw reads, three reference genomes", size=14, italic=True, color=ICE)
text(s, 0.8, 6.45, 11.5, 0.4, "Yusuf Maulana  ·  MBZUAI  ·  draft for ASHG / ISMB 2027  ·  status 25 Oct 2026", size=14, color=ICE)
notes(s, "One line: we take the same raw RNA-seq reads and process them on three reference genomes, and ask what changes "
         "in the genetic map of gene expression.")

# 2 ---------------------------------------------------------------------------------------------------------
s = slide(); title(s, "Why the reference genome matters", "Every step of an eQTL study is anchored to the reference")
stat(s, 0.6, 1.9, 3.6, "151 Mb", "of GRCh38 is unknown sequence (gaps, centromeres, duplications)", color=G38)
stat(s, 0.6, 3.6, 3.6, "238 Mb", "added or corrected by the complete T2T-CHM13 genome", color=T2T)
stat(s, 0.6, 5.3, 3.6, "~20", "known GRCh38 errors patched by FixItFelix", color=FFX)
card(s, 4.7, 1.9, 8.0, 4.85)
text(s, 5.0, 2.1, 7.4, 0.5, "Reads from missing sequence do not vanish", size=20, bold=True, color=NAVY, font=HEAD)
bullets(s, 5.0, 2.8, 7.4, 3.8, [
    ("An eQTL ", "is a DNA variant whose alleles go with higher or lower expression of a nearby gene"),
    ("Genotypes, read placement and gene models ", "all depend on the reference"),
    ("False duplications ", "split one gene's reads between two copies"),
    ("Missing or collapsed copies ", "pile a paralog's reads onto the gene that is present"),
    ("Result: ", "eQTLs that are artefacts, and real ones that are missed"),
], size=17, after=10)
footer(s)
notes(s, "GRCh38 has about 151 Mb of unknown sequence, false duplications such as part of chr21, and collapsed duplications. "
         "T2T-CHM13 adds or corrects 238 Mb (Nurk et al. 2022). When a stretch of genome is missing, reads from it map "
         "to the most similar sequence that exists, and that can create or hide an eQTL.")

# 3 ---------------------------------------------------------------------------------------------------------
s = slide(); title(s, "What is already known, and the gap", "FixItFelix showed artefactual eQTLs, but only inside its patches")
rows = [["", "FixItFelix (Behera et al. 2023)", "This study"],
        ["Regions", "~20 patched GRCh38 regions (a few Mb)", "Whole genome on a complete reference"],
        ["Genotypes", "GRCh38 only", "Native T2T and native GRCh38 panels"],
        ["Gene models", "GRCh38", "RefSeq RS_2025_08 on each assembly"],
        ["Cohort", "Geuvadis, 449 LCLs", "MAGE, 731 LCLs, 26 populations (+ Geuvadis)"],
        ["Key finding", "PWP2, GATD3A eQTLs are artefacts", "57 genes change beyond the patch"]]
table(s, 0.6, 1.9, 12.1, rows, [2.2, 4.6, 5.3], size=15, row_h=0.62, name="Comparison table")
text(s, 0.6, 5.85, 12.1, 0.9, [("Gap: ", {"bold": True, "color": NAVY}),
     ("no eQTL study has been run natively on T2T-CHM13 (literature checked Oct 2026).", {})], size=17)
footer(s)
notes(s, "FixItFelix masked false copies and added decoys for collapsed genes in GRCh38, then remapped 449 Geuvadis samples. "
         "Some GRCh38 eQTLs disappeared, for example PWP2 and GATD3A on chr21. Their analysis covered only the patched "
         "regions and kept GRCh38 genotypes and gene models. We do the whole genome with everything native to T2T.")

# 4 ---------------------------------------------------------------------------------------------------------
s = slide(); title(s, "The question")
card(s, 0.6, 1.45, 12.1, 1.75, fill=NAVY)
text(s, 0.95, 1.6, 11.4, 1.5, [("H0  ", {"bold": True, "color": FFX, "size": 26}),
     ("Most reference-driven eQTL changes lie outside the regions FixItFelix patched", {"color": WHITE, "size": 24})],
     font=HEAD, anchor=MSO_ANCHOR.MIDDLE)
hyp = [("H1", "New variants", "Variants in T2T-only sequence act as eQTLs", NEW),
       ("H2", "New or fixed genes", "Genes absent or mis-assembled in GRCh38 gain eQTLs", T2T),
       ("H3", "Read placement", "Better read assignment changes eQTLs in duplicated genes", SD),
       ("H4", "Artefact removal", "Some GRCh38 eQTLs disappear on T2T", G38)]
for i, (k, h, d, c) in enumerate(hyp):
    x = 0.6 + i * 3.1
    card(s, x, 3.55, 2.85, 2.35)
    dot(s, x + 0.25, 3.8, 0.6, c, k, size=15)
    text(s, x + 0.2, 4.5, 2.5, 0.45, h, size=16, bold=True, color=NAVY)
    text(s, x + 0.25, 4.95, 2.4, 0.9, d, size=14, color=INK)
text(s, 0.6, 6.15, 12.1, 0.7, [("Design: ", {"bold": True, "color": NAVY}),
     ("the patched GRCh38 is a third reference, so what known patches fix is separated from what a complete genome adds.", {})],
     size=16)
footer(s)
notes(s, "The central claim is H0. H1 to H4 are the mechanisms by which a reference can change an eQTL. Using the "
         "FixItFelix-patched GRCh38 as a third reference lets us subtract what is already known.")

# 5 ---------------------------------------------------------------------------------------------------------
s = slide(); title(s, "Data: two lymphoblastoid cohorts", "Same cell type, 1000 Genomes individuals with genotypes on both references")
picture(s, fig(2), 0.5, 1.6, 8.3, 5.35, align="left")
stat(s, 9.2, 1.75, 3.6, "731", "MAGE individuals, 26 populations, 2x150 bp, one lab (discovery)", color=T2T)
stat(s, 9.2, 3.55, 3.6, "451", "Geuvadis individuals, 5 populations, 2x76 bp (replication)", color="56B4E9")
stat(s, 9.2, 5.35, 3.6, "93", "people in both cohorts", color=MUTED)
footer(s)
notes(s, "MAGE (Taylor et al. 2024) has 779 runs from 731 individuals across all 26 1000 Genomes populations, sequenced "
         "in one lab with 150 bp paired reads, about 28 million read pairs per run. Geuvadis has 451 individuals with T2T "
         "genotypes, four European populations plus Yoruba, 76 bp reads. Longer reads help in duplicated regions, which is "
         "why MAGE is the discovery cohort.")

# 6 ---------------------------------------------------------------------------------------------------------
s = slide(); title(s, "Where the data come from", "Public data only; every download verified by checksum or read count")
rows = [["Data", "Source / accession", "Publication", "Used for"],
        ["MAGE RNA-seq", "NCBI SRA PRJNA851328", "Taylor 2024, Nature", "Discovery expression"],
        ["Geuvadis RNA-seq", "ENA PRJEB3366", "Lappalainen 2013, Nature", "Replication expression"],
        ["1KGP genotypes, T2T", "Zenodo 7612953 (SHAPEIT5 v1.1)", "Aganezov 2022; Rhie 2023", "T2T-native genotypes"],
        ["1KGP genotypes, GRCh38", "NYGC 2022 phased panel", "Byrska-Bishop 2022, Cell", "GRCh38 genotypes"],
        ["T2T-CHM13v2.0 + RefSeq", "NCBI GCF_009914755.1", "Nurk 2022, Science", "Reference 1"],
        ["GRCh38 no-alt + RefSeq", "NCBI GCA_000001405.15", "Genome Reference Consortium", "Reference 2"],
        ["FixItFelix patched GRCh38", "Zenodo 7532965", "Behera 2023, Genome Biol", "Reference 3, control"],
        ["T2T-only sequence, SDs", "T2T consortium; UCSC chains", "Nurk 2022; Vollger 2022", "Gene region classes"],
        ["ctyper copy numbers", "Zenodo 13683015", "Ma & Chaisson 2025, Nat Genet", "Planned covariate"]]
table(s, 0.6, 1.75, 12.1, rows, [3.0, 3.6, 3.2, 2.3], size=13, row_h=0.5, name="Data sources table")
text(s, 0.6, 6.85 - 0.3, 12.1, 0.35, "Full list, download commands and checks: workbook sheet 'Data sources'; docs/t2t-eqtl.md section 4",
     size=11, italic=True, color=MUTED)
footer(s)
notes(s, "Everything is public. Reads come from NCBI SRA and are checked with vdb-validate and against ENA's read counts. "
         "The T2T genotypes are the phased 1000 Genomes panel called on CHM13; the GRCh38 genotypes are the NYGC panel. "
         "Gene models are the same NCBI RefSeq release on both assemblies, so genes match by GeneID.")

# 7 ---------------------------------------------------------------------------------------------------------
s = slide(); title(s, "Pipeline: same reads, three references")
picture(s, fig(1), 0.4, 1.15, 7.3, 5.85, align="left")
steps = [("1", "Download and verify", "SRA reads, checksums, read counts"),
         ("2", "Quantify three times", "STAR and salmon per reference (Nextflow)"),
         ("3", "Genotypes", "Native panels; 12.6 M variants matched"),
         ("4", "Map eQTLs", "tensorQTL, 10,000 permutations"),
         ("5", "Compare", "Same gene, same variant; STAR and salmon")]
for i, (k, h, d) in enumerate(steps):
    y = 1.45 + i * 1.07
    dot(s, 7.95, y, 0.5, NAVY, k, size=15)
    text(s, 8.6, y - 0.07, 4.2, 0.42, h, size=17, bold=True, color=NAVY)
    text(s, 8.6, y + 0.33, 4.2, 0.6, d, size=14, color=MUTED)
footer(s)
notes(s, "Each run is downloaded once and aligned three times, so the only thing that changes is the reference. "
         "Expression is quantified with STAR gene counts and, as a check, salmon. Covariates: 5 genotype PCs, sex and "
         "60 expression PCs. The comparison uses variants matched between the two genotype panels.")

# 8 ---------------------------------------------------------------------------------------------------------
s = slide(); title(s, "T2T places reads less ambiguously", "Same 779 MAGE runs aligned to each reference")
picture(s, fig(3), 0.5, 1.7, 8.7, 4.6, align="left")
stat(s, 9.35, 2.0, 3.4, "4.2% → 2.1%", "multi-mapped reads, GRCh38 → T2T (median)", color=T2T, big_size=30)
stat(s, 9.35, 3.9, 3.4, "+0.4 pts", "uniquely mapped reads on T2T (84.9% vs 84.5%)", color=T2T, big_size=30)
footer(s)
notes(s, "Grey lines join the same run on the three references. Halving the multi-mapped reads is the first sign that "
         "read placement changes. The patched GRCh38 behaves almost exactly like GRCh38 here, because the patch is small.")

# 9 ---------------------------------------------------------------------------------------------------------
s = slide(); title(s, "A fair comparison of genotypes", "Only variants present on both references enter the main test")
picture(s, fig(4), 0.6, 1.65, 12.1, 4.55)
card(s, 0.6, 6.25, 12.1, 0.7)
text(s, 0.85, 6.3, 11.6, 0.6, [("12.6 M variants matched ", {"bold": True, "color": NAVY}),
     ("(same alleles, frequency and genotypes in the same people), so any change comes from expression, not genotypes", {})],
     size=15, anchor=MSO_ANCHOR.MIDDLE)
footer(s)
notes(s, "Variants are called independently on each reference. LiftOver is used only to pair them; pairs must have the "
         "same alleles and allele frequency within 0.02, and we check that person-level genotypes agree. 0.13 M pairs "
         "with different frequencies are paralogous positions and are excluded.")

# 10 --------------------------------------------------------------------------------------------------------
s = slide(); title(s, "The overall map is stable", "About 15,000 eGenes on every reference")
picture(s, fig(5), 0.5, 1.7, 8.4, 5.0, align="left")
stat(s, 9.3, 2.0, 3.5, "~15,000", "eGenes (FDR 5%) on each of the three references", color=T2T)
stat(s, 9.3, 3.85, 3.5, "≈ MAGE paper", "which reports ~15,000 eGenes on GRCh38", color=MUTED, big_size=28)
footer(s)
notes(s, "An eGene is a gene with at least one significant eQTL. Totals barely move, and they agree with the MAGE paper, "
         "so the pipeline behaves as expected. The question is which genes change, not how many.")

# 11 --------------------------------------------------------------------------------------------------------
s = slide(); title(s, "Positive control: the known artefacts reappear", "FixItFelix's chr21 genes, mapped in MAGE")
picture(s, fig(6), 0.5, 1.65, 7.4, 5.2, align="left")
pc = pd.read_csv(DD / "positive-control.tsv", sep="\t").pivot(index="symbol", columns="reference", values="qval")
pc = pc.sort_values("grch38")
rows = [["Gene", "T2T", "GRCh38", "Patched"]] + [[g, f"{r['chm13v2.0']:.0e}", f"{r['grch38']:.0e}", f"{r['grch38-ffx']:.0e}"]
                                                  for g, r in pc.iterrows()]
table(s, 8.2, 1.8, 4.5, rows, [1.2, 1.1, 1.1, 1.1], size=13, row_h=0.45, name="Control q-values")
text(s, 8.2, 4.65, 4.5, 0.4, "q-values; smaller = stronger", size=11, italic=True, color=MUTED)
text(s, 8.2, 5.15, 4.5, 1.6, [("The patch changes 9 genes, all in its regions. ", {"bold": True, "color": NAVY}),
     ("T2T agrees with the patched reference.", {})], size=16)
footer(s)
notes(s, "PWP2 has q = 1e-210 on unpatched GRCh38 but about 1e-33 on T2T and on the patched GRCh38. GATD3, CBS, ICOSLG and "
         "U2AF1 behave the same way. This reproduces FixItFelix independently, in a larger and more diverse cohort.")

# 12 --------------------------------------------------------------------------------------------------------
s = slide(); title(s, "Beyond the patch: 57 genes change", "Same gene, same variant: T2T vs patched GRCh38")
picture(s, fig(7), 0.5, 1.6, 6.3, 5.35, align="left")
stat(s, 7.3, 1.7, 2.6, "57", "genes change robustly", color=SD, big_size=54)
stat(s, 10.1, 1.7, 2.6, "0", "of them in FixItFelix regions", color=NAVY, big_size=54)
card(s, 7.3, 3.65, 5.4, 2.6)
text(s, 7.55, 3.8, 5.0, 0.45, "A change counts only if", size=17, bold=True, color=NAVY)
bullets(s, 7.55, 4.3, 5.0, 1.9, ["STAR and salmon both call it, same direction",
                                  "the variant's genotypes agree between panels",
                                  "no annotation, MHC or PAR flag"], size=15, after=6)
text(s, 7.3, 6.4, 5.4, 0.5, "About 7 times as many genes as the patch fixes", size=15, italic=True, color=MUTED)
footer(s)
notes(s, "Each point is a gene: the effect of the same variant on the patched GRCh38 (x) and on T2T (y). Almost all ~15,000 "
         "genes sit on the diagonal. 57 change robustly, none in FixItFelix's regions. An earlier version had 69; 12 were "
         "removed after we found an allele-orientation error for some indels and fixed it using the genotypes themselves.")

# 13 --------------------------------------------------------------------------------------------------------
s = slide(); title(s, "What the 57 genes are", "Gene families and highly expressed genes with pseudogene copies")
picture(s, fig(8), 0.5, 1.6, 6.6, 5.35, align="left")
for i, (big, lab, c) in enumerate([("33", "in segmental duplications (≥ 98% identity)", SD),
                                    ("31", "with opposite-sign effects at the same variant", G38),
                                    ("33", "with one clear gene that takes their reads", T2T)]):
    stat(s, 7.6, 1.65 + i * 1.45, 5.1, big, lab, color=c, big_size=40)
card(s, 7.6, 6.0, 5.1, 0.9)
text(s, 7.8, 6.05, 4.8, 0.8, "CCZ1, ASAH2, PMS2, NBPF10, NPIPB13, TBC1D3E, STAG3L1, SMN2, RPL9, RPSA, UBC",
     size=13, italic=True, color=INK, anchor=MSO_ANCHOR.MIDDLE)
footer(s)
notes(s, "Changes go both ways: 30 genes are stronger on T2T, 27 on the patched GRCh38. About half are in segmental "
         "duplications; the rest include very highly expressed genes such as ribosomal-protein genes, which have many "
         "processed pseudogene copies.")

# 14 --------------------------------------------------------------------------------------------------------
s = slide(); title(s, "Mechanism: reads move between copies", "Per person, one gene loses what its partner gains")
picture(s, fig(9), 0.5, 1.6, 6.9, 5.35, align="left")
mech = [("Paralog swaps", "CCZ1 ↔ CCZ1B (r = −0.96), ASAH2 ↔ ASAH2B, PMS2 ↔ PMS2CL", SD),
        ("Processed pseudogenes", "RPSA → RPSA2 (r = −0.94), H3-3A → H3P6", G38),
        ("Genes only in T2T", "PDPR, PDPR2P → a gene in sequence GRCh38 lacks (r = −0.82)", NEW)]
for i, (h, d, c) in enumerate(mech):
    y = 1.7 + i * 1.75
    card(s, 7.7, y, 5.0, 1.55)
    dot(s, 7.95, y + 0.25, 0.35, c)
    text(s, 8.45, y + 0.15, 4.1, 0.45, h, size=17, bold=True, color=NAVY)
    text(s, 8.45, y + 0.6, 4.1, 0.9, d, size=14, color=INK)
footer(s)
notes(s, "For each gene we take, per person, the change in expression between the references and find the gene whose "
         "change mirrors it. CCZ1 and CCZ1B swap reads almost exactly. RPSA loses reads to its processed pseudogene RPSA2. "
         "PDPR on GRCh38 holds part of the reads, about 8% of its change, that T2T assigns to a gene in sequence absent from "
         "GRCh38; for PDPR2P it is 95%. Missing sequence creating "
         "an eQTL artefact outside every patch.")

# 15 --------------------------------------------------------------------------------------------------------
s = slide(); title(s, "Same people, same variant, different answer", "Expression by genotype on two references")
picture(s, fig(10), 0.6, 1.65, 12.1, 4.75)
cols = [("PWP2", "the known artefact: strong only on unpatched GRCh38"), ("CCZ1", "the effect reverses direction"),
        ("SMN2", "much stronger on T2T")]
for i, (g, d) in enumerate(cols):
    text(s, 0.6 + i * 4.1, 6.4, 3.9, 0.55, [(g + ": ", {"bold": True, "color": NAVY}), (d, {})], size=15)
footer(s)
notes(s, "The reference can create, remove or reverse an eQTL. These genes are tested in GTEx and used in colocalisation "
         "and fine-mapping, so the choice of reference changes downstream conclusions.")

# 16 --------------------------------------------------------------------------------------------------------
s = slide(); title(s, "Limits and next steps")
for x, head, items, c in [(0.6, "Caveats", [
        "salmon alone calls ~550 genes changed (transcript sets differ between assemblies), so both methods are required",
        "An indel allele-orientation error added 12 false genes; fixed by taking orientation from the genotypes",
        "Unmapped 'other' reads: 0.93% on T2T vs 0.09%, not yet explained",
        "24 of 57 genes have no single read partner yet"], G38),
    (6.8, "Next", [
        "Replication in Geuvadis (running, about 27 Oct)",
        "eQTLs from T2T-only variants (running)",
        "chrX with reference-specific PAR boundaries",
        "Paralog copy number (ctyper) as a covariate",
        "GM12878 regulatory elements and GWAS colocalisation"], T2T)]:
    card(s, x, 1.45, 5.9, 5.4)
    dot(s, x + 0.3, 1.7, 0.4, c)
    text(s, x + 0.85, 1.65, 4.8, 0.5, head, size=22, bold=True, color=NAVY, font=HEAD)
    bullets(s, x + 0.3, 2.35, 5.35, 4.4, items, size=15, after=9)
footer(s)
notes(s, "Be upfront about the caveats. Salmon distributes reads across transcript models, and RefSeq's transcript sets "
         "differ between assemblies, so salmon alone picks up annotation differences; requiring both quantifiers removes "
         "them. Replication in Geuvadis is the main remaining test.")

# 17 --------------------------------------------------------------------------------------------------------
s = slide(dark=True); title(s, "Conclusions", dark=True)
concl = [("1", "The overall eQTL map is stable across references (~15,000 eGenes each)"),
         ("2", "Known GRCh38 artefacts reappear and T2T fixes them, as the patch does"),
         ("3", "Beyond the patch, T2T changes 57 genes, mostly by moving reads between paralogs, pseudogenes and new genes")]
for i, (k, t) in enumerate(concl):
    y = 1.7 + i * 1.45
    dot(s, 0.8, y, 0.7, [T2T, FFX, SD][i], k, size=20)
    text(s, 1.8, y - 0.05, 10.8, 0.85, t, size=22, color=WHITE, anchor=MSO_ANCHOR.MIDDLE)
text(s, 0.8, 6.15, 11.8, 0.6, "Take-home: check gene families on a complete reference before trusting their eQTLs",
     size=18, italic=True, color=ICE)
notes(s, "Close on the take-home. Mention that code and results will be released and that Geuvadis replication is under way.")

# 18 --------------------------------------------------------------------------------------------------------
s = slide(); title(s, "Backup: key numbers", "MAGE, 731 individuals, chr1-22, matched variants")
eg = pd.read_csv(DD / "egenes.tsv", sep="\t").pivot(index="reference", columns="method", values="egenes")
rows = [["Reference", "eGenes, STAR", "eGenes, salmon"]] + [[n, f"{eg.loc[r, 'star']:,}", f"{eg.loc[r, 'salmon']:,}"]
        for r, n in [("chm13v2.0", "T2T-CHM13"), ("grch38", "GRCh38"), ("grch38-ffx", "Patched GRCh38")]]
table(s, 0.6, 1.75, 5.6, rows, [2.2, 1.7, 1.7], size=13, row_h=0.55, name="eGene table")
rows = [["Genes changed", "STAR", "salmon", "Both", "Clean"],
        ["Patched vs GRCh38", 9, 9, 8, 6], ["T2T vs GRCh38", 145, 563, 79, 63], ["T2T vs patched", 138, 549, 71, 57]]
table(s, 6.6, 1.75, 6.1, rows, [2.0, 1.05, 1.05, 1.0, 1.0], size=13, row_h=0.55, name="Change table")
text(s, 0.6, 4.2, 12.1, 0.4, "Clean = no annotation, MHC or PAR flag. Bonferroni 5% on the same-variant effect difference.",
     size=12, italic=True, color=MUTED)
for i, (big, lab) in enumerate([("779", "MAGE runs, each aligned to 3 references"), ("6,157", "CPU hours for MAGE expression"),
                                  ("12.6 M", "variants matched between panels"), ("~44", "CPU hours per eQTL scan (reference x method)")]):
    stat(s, 0.6 + i * 3.1, 4.85, 2.8, big, lab, color=NAVY, big_size=34)
footer(s)
notes(s, "Numbers for questions. Full tables are in the workbook (Key results, 57 genes).")

prs.core_properties.title = "A complete reference genome changes the eQTL map beyond known GRCh38 fixes"
prs.core_properties.author = "Yusuf Maulana"
OUT.parent.mkdir(parents=True, exist_ok=True)
prs.save(OUT)
print(OUT, state["n"], "slides")
