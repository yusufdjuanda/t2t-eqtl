"""Presentation: the oral-talk deck (PowerPoint) built from the figures and the workbook's content.

    pixi run python scripts/07-presentation/build-deck.py

Inputs: results/07-presentation/figures/*.png (pipeline.dot, make-figures.R) and results/07-presentation/data/*.tsv.
Output: results/07-presentation/t2t-eqtl-talk-v2.pptx (16:9, 25 slides + 13 backup; speaker notes on every slide,
in three parts: what the slide shows, why it matters, where it fits in the story).
Story (v2, sections set 2026-11-03): introduction -> data and design -> 1 expression -> 2 genotypes ->
3 eQTL mapping -> 4 comparison (methods; how many change, with the patch as sanity check; top-10 table; segmental
duplications; reads moving) -> 5 link to GWAS (one slide) -> limits, conclusions (simplified 2026-11-03; replication,
copy number, examples and allele-specific expression are backups). Each section opens with a methods slide (method_slide()), then its results. A section label
at the top right of each slide shows where it is in the story.
The 2026-10-29 deck (positive control first) is t2t-eqtl-talk.pptx, built by this script's earlier version in git.
Text follows build-workbook.py (Talk outline) and docs/t2t-eqtl.md; update them together.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm
from PIL import Image as PILImage
from pptx import Presentation
from pptx.chart.data import CategoryChartData
from pptx.enum.chart import XL_CHART_TYPE, XL_LEGEND_POSITION
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.oxml.ns import qn
from pptx.oxml.xmlchemy import OxmlElement
from pptx.util import Emu, Inches, Pt

ROOT = Path(__file__).resolve().parents[2]
FD = ROOT / "results/07-presentation/figures"
DD = ROOT / "results/07-presentation/data"
OUT = ROOT / "results/07-presentation/t2t-eqtl-talk-v2.pptx"

HEAD, BODY = "Cambria", "Calibri"
NAVY, INK, MUTED, TINT, WHITE = "14213D", "1F2933", "5B6770", "EEF2F7", "FFFFFF"
T2T, G38, FFX, SD, NEW = "0072B2", "D55E00", "E69F00", "CC79A7", "009E73"
ICE = "C9D6EA"
rgb = lambda h: RGBColor.from_string(h)


def reversed_again(d):
    """Slopes again of opposite sign under STAR and salmon, none near zero (|slope| >= 0.05; TCAF2 and ELOB are one-sided)."""
    opp = lambda m: (np.sign(d[f"{m}_slope_a"]) != np.sign(d[f"{m}_slope_b"])) & (d[f"{m}_slope_a"].abs() >= 0.05) & (d[f"{m}_slope_b"].abs() >= 0.05)
    return opp("star") & opp("salmon")

prs = Presentation()
prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
BLANK = prs.slide_layouts[6]
W, H = 13.333, 7.5
state = {"n": 0}


def slide(dark=False):
    s = prs.slides.add_slide(BLANK)
    s.background.fill.solid(); s.background.fill.fore_color.rgb = rgb(WHITE)   # white everywhere (decided 2026-11-02)
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
    if state.get("sec") and not dark:   # where we are in the story (section label, top right)
        text(s, 8.7, 0.05, 4.0, 0.28, state["sec"].upper(), size=11, bold=True, color=MUTED, align=PP_ALIGN.RIGHT,
             name="Section")
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


def gt_table(s, x, y, w, rows, col_w, size=12, row_h=0.36, align=None, bold_cols=(0,), colors=None, name=None):
    """A clean booktabs-style table (as R gt): white, no vertical lines, rules above and below the header and at
    the bottom. align: one of 'l', 'r', 'c' per column. colors: {(row, col): hex} overrides for text colour."""
    t = s.shapes.add_table(len(rows), len(rows[0]), Inches(x), Inches(y), Inches(w), Inches(row_h * len(rows)))
    if name:
        t.name = name
    t = t.table
    t.first_row = False; t.horz_banding = False
    al = {"l": PP_ALIGN.LEFT, "r": PP_ALIGN.RIGHT, "c": PP_ALIGN.CENTER}
    align = align or ["l"] * len(rows[0])
    colors = colors or {}

    def ln(tag, w_pt, color=INK):
        e = OxmlElement(tag)
        if w_pt:
            e.set("w", str(int(w_pt * 12700)))
            f = OxmlElement("a:solidFill"); c = OxmlElement("a:srgbClr"); c.set("val", color); f.append(c); e.append(f)
        else:
            e.append(OxmlElement("a:noFill"))
        return e

    for j, cw in enumerate(col_w):
        t.columns[j].width = Inches(cw)
    last = len(rows) - 1
    for i, r in enumerate(rows):
        t.rows[i].height = Inches(row_h)
        for j, v in enumerate(r):
            c = t.cell(i, j); c.text = ""
            c.margin_left = c.margin_right = Inches(0.07); c.margin_top = c.margin_bottom = Inches(0.02)
            c.vertical_anchor = MSO_ANCHOR.MIDDLE
            c.fill.solid(); c.fill.fore_color.rgb = rgb(WHITE)
            tcPr = c._tc.get_or_add_tcPr()
            for k in list(tcPr):
                if k.tag in (qn("a:lnL"), qn("a:lnR"), qn("a:lnT"), qn("a:lnB")):
                    tcPr.remove(k)
            top = 1.5 if i == 0 else (0.75 if i == 1 else 0)
            bot = 1.5 if i == last else (0.75 if i == 0 else 0.25)
            bc = INK if (i == 0 or i == last) else "D5DAE0"
            for n, e in enumerate([ln("a:lnL", 0), ln("a:lnR", 0), ln("a:lnT", top), ln("a:lnB", bot, bc)]):
                tcPr.insert(n, e)
            pp = c.text_frame.paragraphs[0]; pp.alignment = al[align[j]]
            run(pp, str(v), size, MUTED if i == 0 else colors.get((i, j), INK), bold=(i == 0 or j in bold_cols),
                italic=(i > 0 and j == 0))
    return t


def notes(s, t):
    s.notes_slide.notes_text_frame.text = t


def notes3(s, what, why, story):
    """Speaker notes in three parts: what the slide shows, why it matters, and where it fits in the story."""
    notes(s, f"WHAT THIS SLIDE SHOWS\n{what}\n\nWHY IT MATTERS\n{why}\n\nWHERE IT FITS IN THE STORY\n{story}")


def method_slide(s, steps, tools, output, span=5.2):
    """Methods slide: numbered steps on the left; tools (with versions) and the output on the right."""
    for i, (h, d) in enumerate(steps):
        y = 1.65 + i * (span / len(steps))
        dot(s, 0.6, y, 0.48, NAVY, str(i + 1), size=15)
        text(s, 1.25, y - 0.06, 6.6, 0.42, h, size=17, bold=True, color=NAVY)
        text(s, 1.25, y + 0.36, 6.6, 0.75, d, size=13.5, color=INK)
    card(s, 8.3, 1.6, 4.4, 3.05)
    text(s, 8.55, 1.72, 4.0, 0.4, "Tools", size=16, bold=True, color=NAVY)
    bullets(s, 8.55, 2.15, 4.0, 2.45, tools, size=13, after=4)
    card(s, 8.3, 4.85, 4.4, 2.0)
    text(s, 8.55, 4.97, 4.0, 0.4, "Output", size=16, bold=True, color=NAVY)
    text(s, 8.55, 5.4, 4.0, 1.4, output, size=13, color=INK)


def fig_bullets(s, png, items, fig_w=7.2, top=1.6, h=5.25, size=15, after=12):
    """One figure on the left, short bullets (bold lead, plain rest) on the right."""
    px, py, pw, ph = picture(s, png, 0.5, top, fig_w, h, align="left")
    x = 0.5 + pw + 0.4
    bullets(s, x, top + 0.05, 12.75 - x, h, items, size=size, after=after, name="Explanation")
    return px, py, pw, ph


fig = lambda n: next(FD.glob(f"fig{n:02d}-*.png")) if isinstance(n, int) else FD / f"{n}.png"

GS = ROOT / "results/06-compare/eqtl-comparison/mage/gene-sets"
PAIRS = [("chm13v2.0-vs-grch38", "T2T vs GRCh38"), ("grch38-ffx-vs-grch38", "Patched vs GRCh38"),
         ("chm13v2.0-vs-grch38-ffx", "T2T vs patched GRCh38")]
sets = {p: pd.read_csv(GS / f"changed-{p}.tsv", sep="\t") for p, _ in PAIRS}
main = {p: t[t.tier.isin(["core", "main"])] for p, t in sets.items()}
in_ffx = {p: int((t.region_class == "ffx_patched").sum()) for p, t in main.items()}
tg, pg, tp = (main[p] for p, _ in PAIRS)
tg_out = set(tg.loc[tg.region_class != "ffx_patched", "symbol"])
kept = len(tg_out & set(tp.symbol))
tg_rep = tg[tg.geuvadis_same_sign.notna()]
tp_rep = tp[tp.geuvadis_same_sign.notna()]


def bar_chart(s, x, y, w, h, cats, series, colors, name, stacked=True, legend=True, size=14):
    """Native (editable) horizontal bar chart; series = [(name, values)]."""
    cd = CategoryChartData(); cd.categories = cats
    for n, v in series:
        cd.add_series(n, v)
    kind = XL_CHART_TYPE.BAR_STACKED if stacked else XL_CHART_TYPE.BAR_CLUSTERED
    gf = s.shapes.add_chart(kind, Inches(x), Inches(y), Inches(w), Inches(h), cd); gf.name = name
    ch = gf.chart
    ch.font.size = Pt(size); ch.font.name = BODY; ch.font.color.rgb = rgb(INK)
    ch.has_legend = legend
    if legend:
        ch.legend.position = XL_LEGEND_POSITION.BOTTOM; ch.legend.include_in_layout = False
    ca = ch.category_axis; ca.reverse_order = True; ca.format.line.color.rgb = rgb(ICE)
    ca.has_major_gridlines = False
    va = ch.value_axis; va.has_major_gridlines = True; va.major_gridlines.format.line.color.rgb = rgb(TINT)
    va.format.line.fill.background(); va.tick_labels.font.color.rgb = rgb(MUTED)
    pl = ch.plots[0]; pl.gap_width = 60; pl.overlap = 100 if stacked else 0
    for ser, c in zip(pl.series, colors):
        ser.format.fill.solid(); ser.format.fill.fore_color.rgb = rgb(c)
        ser.data_labels.show_value = True; ser.data_labels.font.size = Pt(size); ser.data_labels.font.bold = True
        ser.data_labels.font.color.rgb = rgb(WHITE)
    return ch


# 1 -------------------------------------------------------------------------------------------------------------
state['sec'] = None
s = slide()
text(s, 0.8, 1.4, 11.5, 2.3, "A complete reference genome changes the eQTL map beyond known GRCh38 fixes",
     size=40, bold=True, font=HEAD, color=NAVY, name="Title")
text(s, 0.8, 3.85, 11.5, 0.6, "T2T-native cis-eQTL mapping in 731 MAGE lymphoblastoid cell lines", size=20, color=MUTED)
for i, (c, lab) in enumerate([(T2T, "T2T-CHM13"), (G38, "GRCh38"), (FFX, "Patched GRCh38")]):
    dot(s, 0.8 + i * 2.6, 4.85, 0.32, c)
    text(s, 1.2 + i * 2.6, 4.82, 2.2, 0.4, lab, size=15, color=INK)
text(s, 0.8, 5.35, 11, 0.4, "Same raw reads, three reference genomes", size=14, italic=True, color=MUTED)
text(s, 0.8, 6.45, 11.5, 0.4, "Yusuf Maulana  ·  MBZUAI  ·  draft for ASHG / ISMB 2027  ·  status 2 Nov 2026", size=14, color=MUTED)
notes3(s,
    "The title and the setup in one line: we took the same raw RNA-seq reads from 731 people and processed them three "
    "times, on GRCh38, on GRCh38 with its known errors patched, and on the complete T2T-CHM13 genome, and asked how the "
    "genetic map of gene expression changes.",
    "eQTL maps are how we read GWAS: they tell us which gene a disease variant probably acts through. Almost every eQTL "
    "map, GTEx, Geuvadis and MAGE included, was built on GRCh38. If the reference itself shapes the map, every study "
    "that uses these catalogues inherits it.",
    "The talk follows the data in the order it was processed: expression, genotypes, eQTLs, then the comparison between "
    "references with checks that the changes are real, and finally the link to GWAS. Each section starts with a short "
    "methods slide, then the results.")

# 2 -------------------------------------------------------------------------------------------------------------
state['sec'] = 'Introduction'
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
notes3(s,
    "What an eQTL needs and why each part depends on the reference. An eQTL is a DNA variant whose alleles go with more "
    "or less RNA from a nearby gene. To find one you need genotypes, expression and gene models, and all three are "
    "defined on a reference genome. GRCh38 has about 151 Mb of unknown sequence, false duplications such as part of "
    "chr21, and collapsed duplications. T2T-CHM13 adds or corrects 238 Mb (Nurk et al. 2022).",
    "Reads from sequence that is missing from the reference do not disappear: the aligner puts them on the most similar "
    "sequence that does exist, or throws them away as ambiguous. Either way the expression of a gene can be wrong, and if "
    "how wrong depends on a person's genotype, a fake eQTL appears or a real one is hidden.",
    "This is the problem. The next slide asks whether anyone has measured how big it is for eQTLs.")

# 3 -------------------------------------------------------------------------------------------------------------
state['sec'] = 'Introduction'
s = slide(); title(s, "The gap and the question", "A complete reference exists, but the eQTL map is still drawn on GRCh38")
gap = [("GRCh38 ", "has gaps, false duplications and collapsed copies"),
       ("T2T-CHM13 ", "is complete, and 1000 Genomes genotypes are already called on it"),
       ("Yet eQTL maps ", "(GTEx, Geuvadis, MAGE) are built on GRCh38 or older"),
       ("No study ", "has mapped eQTLs natively on T2T (literature checked Oct 2026)")]
card(s, 0.6, 1.8, 5.9, 4.1)
text(s, 0.9, 1.95, 5.4, 0.5, "The gap", size=20, bold=True, color=NAVY, font=HEAD)
bullets(s, 0.9, 2.55, 5.4, 3.3, gap, size=16, after=10)
card(s, 6.85, 1.8, 5.85, 4.1)
text(s, 7.15, 1.95, 5.3, 0.5, "Research question", size=20, bold=True, color=NAVY, font=HEAD)
text(s, 7.15, 2.55, 5.3, 1.9, "When the same reads are mapped to the complete T2T genome instead of GRCh38, which eQTLs change, "
     "and are the changes more than the known, patchable GRCh38 errors?", size=19, color=INK, font=HEAD)
text(s, 7.15, 4.55, 5.3, 1.2, [("H0  ", {"bold": True, "color": NAVY}),
     ("Most reference-driven eQTL changes lie outside the regions FixItFelix patched", {"color": INK})], size=16)
text(s, 0.6, 6.15, 12.1, 0.7, [("Plan: ", {"bold": True, "color": NAVY}),
     ("compare T2T with GRCh38, check the known patch on its own, then subtract it: T2T vs patched GRCh38.", {})], size=16)
footer(s)
notes3(s,
    "The gap and our question. GRCh38's errors are known and T2T has existed since 2022, with the 1000 Genomes "
    "genotypes already called on it, yet eQTL maps are still built on GRCh38, and our literature searches (Oct 2026) "
    "found no eQTL study run natively on T2T.",
    "The question is not only whether T2T changes eQTLs, but whether it changes them beyond the errors we already know "
    "how to patch. H0 makes that testable: if most changes lie outside the patched regions, patching GRCh38 is not "
    "enough and the field needs a complete reference.",
    "The plan has three comparisons: T2T vs GRCh38, the patch on its own as a check, and T2T vs patched GRCh38. Next: "
    "what earlier work did, which explains why there is a patched GRCh38 at all. The four mechanisms H1 to H4 are in the backup.")

# 4 -------------------------------------------------------------------------------------------------------------
state['sec'] = 'Introduction'
s = slide(); title(s, "Has anyone done this?", "T2T has been used for variants and regulation; for eQTLs only as a patch on GRCh38")
rows = [["Work", "What it did on T2T", "eQTLs on T2T?"],
        ["Aganezov 2022, Science", "1000 Genomes variants called on T2T", "No"],
        ["Gershman 2022, Science", "ENCODE ChIP-seq remapped to T2T", "No"],
        ["Ma & Chaisson 2025, Nat Genet", "Paralog copy number (ctyper)", "No"],
        ["Morrissey 2026, bioRxiv", "Regulation across duplications", "No"],
        ["Behera 2023, Genome Biol", "FixItFelix: T2T as a patch for GRCh38", "Patched regions only"]]
t = table(s, 0.6, 1.8, 6.6, rows, [2.5, 2.75, 1.35], size=13, row_h=0.6, name="Prior work table")
for j in range(3):
    c = t.cell(5, j); c.fill.solid(); c.fill.fore_color.rgb = rgb("FDF0D5")
card(s, 7.55, 1.8, 5.15, 4.15)
text(s, 7.85, 1.95, 4.7, 0.5, "FixItFelix (Behera et al. 2023)", size=19, bold=True, color=NAVY, font=HEAD)
bullets(s, 7.85, 2.55, 4.65, 2.35, [
    ("Uses T2T ", "to fix ~20 known GRCh38 errors: false copies masked, missing copies added"),
    ("Reference stays GRCh38, ", "with GRCh38 genotypes and gene models"),
    ("Geuvadis, 449 LCLs: ", "PWP2 and GATD3A eQTLs are artefacts"),
], size=14, after=7)
for i, (c, lab) in enumerate([(G38, "GRCh38"), (FFX, "+ patch"), (T2T, "T2T")]):
    dot(s, 7.9 + i * 1.6, 5.1, 0.28, c)
    text(s, 8.25 + i * 1.6, 5.05, 1.3, 0.4, lab, size=13, color=INK)
text(s, 7.85, 5.45, 4.7, 0.45, "They tested the first step; we test the whole way", size=13, italic=True, color=MUTED)
text(s, 0.6, 6.15, 12.1, 0.7, [("Open: ", {"bold": True, "color": NAVY}),
     ("everything outside the patched regions. We add patched GRCh38 as a third reference to separate the two.", {})], size=16)
footer(s)
notes3(s,
    "Prior work. Variants, ChIP-seq, copy number and duplications have all been studied on T2T, but no eQTL map. The "
    "closest work is FixItFelix (Behera et al. 2023). The errors themselves were found earlier, by Wagner et al. 2022 "
    "(chr21 false duplications) and Aganezov et al. 2022 (more false duplications and collapses, by comparing GRCh38 "
    "with T2T). FixItFelix masked the false copies, added the missing copies from T2T as decoys, and in 449 Geuvadis "
    "samples showed that eQTLs such as PWP2 and GATD3A are artefacts.",
    "FixItFelix chose to patch GRCh38 rather than move to T2T, to keep GRCh38 coordinates and resources such as gnomAD "
    "and to fix existing alignments in minutes. That was practical, but it can only fix errors someone has already "
    "found: about 20 regions, a few Mb, with GRCh38 genotypes and gene models.",
    "So we use their patched GRCh38 as our third reference. Whatever still differs between it and T2T is something the "
    "known fix does not cover. Next: the journey of the data through the study.")

# ROADMAP -------------------------------------------------------------------------------------------------------
state['sec'] = 'Introduction'
s = slide(); title(s, "The journey of the data", "Same people, same reads; only the reference genome changes")
road = [("1", "Expression", "Where do the reads go on each reference?", T2T),
        ("2", "Genotypes", "Is it the same variant on both references?", NEW),
        ("3", "eQTL mapping", "Which variants affect which genes?", FFX),
        ("4", "Comparison", "Does the same variant have the same effect? Is it real?", G38),
        ("5", "Link to GWAS", "Does it change the gene behind a disease signal?", SD)]
for i, (k, h, d, c) in enumerate(road):
    x = 0.6 + i * 2.46
    card(s, x, 1.9, 2.25, 3.2)
    dot(s, x + 0.2, 2.1, 0.62, c, k, size=20)
    text(s, x + 0.2, 2.95, 2.0, 0.5, h, size=17, bold=True, color=NAVY)
    text(s, x + 0.2, 3.5, 1.9, 1.5, d, size=14, color=INK)
    if i < 4:
        text(s, x + 2.2, 3.2, 0.3, 0.5, "›", size=28, color=MUTED, align=PP_ALIGN.CENTER)
card(s, 0.6, 5.45, 12.1, 1.3)
text(s, 0.9, 5.5, 11.5, 1.2, [("Design rule: ", {"bold": True, "color": NAVY}),
     ("same individuals, same raw reads, same software and versions, same gene annotation release, same covariates. "
      "Any difference can then only come from the reference.", {})], size=17, anchor=MSO_ANCHOR.MIDDLE)
footer(s)
notes3(s,
    "The roadmap of the talk. It follows the data the way it was processed: (1) RNA-seq reads are aligned and counted on "
    "each reference; (2) genotypes are called on each reference and paired so we test the same variant; (3) eQTLs are "
    "mapped on each reference; (4) the same gene and variant are compared between references, with checks that the changes are real; (5) "
    "we ask whether they change the gene a disease signal points to.",
    "The design rule is what makes the comparison meaningful: same people, same reads, same software, same annotation "
    "release, same covariates. Much of the pipeline exists only to remove differences that are not the reference "
    "(contigs, chrY masking, gene matching by NCBI GeneID, shared genotype PCs, paired variants). If anything else "
    "differed, a change could not be blamed on the reference.",
    "The section label at the top right of each slide shows which step we are in. First, the data.")

# 5 -------------------------------------------------------------------------------------------------------------
state['sec'] = 'Data and design'
s = slide(); title(s, "Data: two lymphoblastoid cohorts", "Same cell type; 1000 Genomes people with genotypes on both references")
fig_bullets(s, fig(2), [
    ("MAGE (discovery): ", "731 people from all 26 1000 Genomes populations; one lab, 2×150 bp reads"),
    ("Geuvadis (replication): ", "451 people, mostly European; several labs, 2×76 bp reads; 93 are also in MAGE"),
    ("Why these cohorts: ", "the same cell type, and genotypes already called on both T2T and GRCh38"),
    ("Longer reads ", "place better in duplications, so MAGE discovers and Geuvadis checks")])
footer(s)
notes3(s,
    "The two cohorts. MAGE (Taylor et al. 2024): 779 RNA-seq runs from 731 people in all 26 1000 Genomes populations, "
    "one lab, 150 bp paired reads, about 28 million read pairs per run. Geuvadis: 451 people with T2T genotypes, four "
    "European populations plus Yoruba, 76 bp reads, several labs. 93 people are in both.",
    "Both are lymphoblastoid cell lines from 1000 Genomes individuals, whose genotypes have been called on both GRCh38 "
    "and T2T. That is what makes a T2T-native eQTL study possible at all. MAGE is the discovery cohort because it is "
    "larger, more diverse and has longer reads, which place better in duplicated regions. Geuvadis is the independent "
    "replication, and it is also the cohort FixItFelix used, so their results can be checked directly.",
    "Next: how the data were processed, with the reference as the only thing that changes.")

# DESIGN --------------------------------------------------------------------------------------------------------
state['sec'] = 'Data and design'
s = slide(); title(s, "Design: one cohort, three references", "Each MAGE run is processed three times from the same reads")
px, py, pw, ph = picture(s, FD / "method-overview.png", 0.6, 1.45, 12.1, 4.35)
text(s, 0.6, 5.9, 12.1, 0.5, [("Patched vs GRCh38 ", {"bold": True, "color": FFX}), ("checks the pipeline (only the patch differs).  ", {}),
     ("T2T vs patched GRCh38 ", {"bold": True, "color": T2T}), ("is the main question: what a complete genome adds beyond the known fix.", {})],
     size=14)
text(s, 0.6, 6.4, 12.1, 0.55, [("Two ways of counting: ", {"bold": True, "color": NAVY}),
     ("STAR counts only reads placed in one location; salmon shares ambiguous reads between transcripts. They fail "
      "differently, so a change must appear with both.", {})], size=14)
footer(s)
notes3(s,
    "The method in one picture. The same 731 people provide RNA-seq and genotypes. The RNA is aligned and counted with "
    "STAR and salmon, and eQTLs are mapped with tensorQTL (66 covariates: 5 genotype PCs, sex, 60 expression PCs), "
    "repeated on each of the three references. Genotypes are filtered with plink2 and variants are matched between "
    "references with liftOver plus allele checks. Each pair of references is then compared at the same gene and variant.",
    "The three references give three pairs with different jobs. GRCh38 vs patched GRCh38 differ only inside the patch, so "
    "this pair must reproduce FixItFelix: it tests the pipeline. T2T vs patched GRCh38 is the main question: everything "
    "that changes there, the known fix does not explain. A change has to appear with both STAR and salmon, two very "
    "different ways of counting, to rule out quantification artefacts. STAR (the GTEx-style standard) counts a read for a "
    "gene only if it lands in one place and overlaps that gene alone; reads that fit several places are discarded, so a "
    "gene can lose reads just because the reference has another copy. salmon keeps those reads and shares them between "
    "matching transcripts by likelihood, so it loses fewer reads but depends on the transcript models, which differ "
    "somewhat between assemblies. A counting artefact shows up in one method only (in the pilot, S100A6 changed 10-fold "
    "in STAR and not at all in salmon). salmon here runs on STAR's alignments, so the two differ in counting, not "
    "alignment.",
    "Now we follow the data, starting with step 1: where the RNA-seq reads go on each reference.")

# M_EXPR --------------------------------------------------------------------------------------------------------
state['sec'] = '1 · Expression'
s = slide(); title(s, "Methods: from reads to gene counts", "The same RNA-seq reads, processed once per reference")
method_slide(s, [
    ("Download and check", "779 MAGE runs from NCBI SRA; checksums (vdb-validate) and read counts against ENA"),
    ("Align to each reference", "STAR, ENCODE settings; T2T-CHM13v2.0, GRCh38 and patched GRCh38; chrM and EBV added, chrY PARs masked"),
    ("Count reads per gene, two ways", "STAR: uniquely placed reads only. salmon: ambiguous reads shared between transcripts"),
    ("Same genes on every reference", "NCBI RefSeq RS_2025_08 on both assemblies, genes matched by GeneID; runs summed per person (731)")],
    [("STAR ", "2.7.11b"), ("salmon ", "1.11.4 (alignment mode)"), ("samtools, SRA Toolkit", ""), ("Nextflow ", "pipeline, Slurm")],
    "Gene counts and TPM for 731 people × 3 references × 2 methods. A gene's expression change counts only if it is ≥ 2-fold with both methods.")
footer(s)
notes3(s,
    "How the expression data were made. Each MAGE run was downloaded once, checked, and aligned three times, once per "
    "reference, with the same software and settings. Reads were counted in two independent ways, STAR and salmon. The "
    "same RefSeq annotation release is used on both assemblies, so a gene is the same gene (same NCBI GeneID) on every "
    "reference. Runs from the same person were summed, giving 731 people.",
    "Everything that is not the reference is held fixed here, so any difference later can be blamed on the reference. "
    "Two counting methods guard against differences that come from how reads are counted rather than where they go.",
    "Next: the results of this step, starting with how the reads were placed.")

# 8 -------------------------------------------------------------------------------------------------------------
state['sec'] = '1 · Expression'
s = slide(); title(s, "T2T places reads less ambiguously", "Same 779 MAGE runs aligned to each reference")
qc = pd.read_csv(DD / "mapping-qc.tsv", sep="\t").median(numeric_only=True)
fig_bullets(s, fig(3), [
    ("Each dot ", "is one MAGE run, aligned to all three references from the same reads"),
    ("Multi-mapped reads halve on T2T: ", f"{qc['multi_pct_grch38']:.1f}% → {qc['multi_pct_chm13v2.0']:.1f}% "
     f"(patched GRCh38 {qc['multi_pct_grch38-ffx']:.1f}%)"),
    ("Uniquely mapped reads ", f"rise slightly ({qc['unique_pct_chm13v2.0']:.1f}% vs {qc['unique_pct_grch38']:.1f}%)"),
    ("Why it matters: ", "gene counts drop multi-mapped reads, so T2T can count reads that GRCh38 throws away")])
footer(s)
notes3(s,
    "Mapping of the 779 MAGE runs. Panel a: the share of reads placed in exactly one location, slightly higher on T2T "
    "(+0.4 points). Panel b: the share placed in several equally good locations, halved on T2T (4.3% to 2.1%), while "
    "the patched GRCh38 (4.2%) is almost the same as GRCh38. Each dot is one run, aligned three times from the same reads.",
    "STAR's gene counts use only uniquely placed reads; multi-mapped reads are thrown away. On GRCh38 these are reads "
    "that land in duplicated or missing sequence, so halving them means T2T can count reads that GRCh38 discards. The "
    "patch fixes its own regions but hardly changes this genome-wide: H0 in miniature, before any eQTL is mapped.",
    "Reads are placed differently. Next: how reads are turned into gene counts, with two different methods.")

# 8a ------------------------------------------------------------------------------------------------------------
state['sec'] = '1 · Expression'
ecs = pd.read_csv(DD / "expression-compare.tsv", sep="\t")
es = ecs[ecs.pair == "chm13v2.0-vs-grch38-ffx"]
n_both = int(es.changed.sum())
n_star = int(((es.log2fc_star.abs() >= 1) & ~es.changed).sum())
n_sal = int(((es.log2fc_salmon.abs() >= 1) & ~es.changed & (es.log2fc_star.abs() < 1)).sum())
s = slide(); title(s, "Two ways of counting reads: STAR and salmon", "Both methods see the same large changes; a change seen by only one is not counted")
picture(s, fig(16), 0.5, 1.6, 7.9, 4.6, align="left")
items = [("STAR ", "counts a read for a gene only if it lands in one place and on that gene alone; ambiguous reads are thrown away", T2T),
         ("salmon ", "keeps ambiguous reads and shares them between matching transcripts by likelihood", NEW),
         ("They fail differently: ", "STAR changes when reads become ambiguous, salmon when transcript models differ", NAVY)]
for i, (h, d, c) in enumerate(items):
    y = 1.7 + i * 1.45
    card(s, 8.7, y, 4.0, 1.3)
    text(s, 8.9, y + 0.1, 3.7, 1.15, [(h, {"bold": True, "color": c}), (d, {})], size=13)
text(s, 0.6, 6.3, 12.1, 0.65, [("T2T vs patched GRCh38: ", {"bold": True, "color": NAVY}),
     (f"{n_both} genes change ≥ 2-fold with both; {n_star} with STAR only and {n_sal} with salmon only are not counted as changes. All genes stay in the eQTL analysis.", {})], size=15)
footer(s)
notes3(s,
    "How reads become gene counts, and why we use two methods. Panel a: each gene's expression change between T2T and "
    f"patched GRCh38 measured with STAR (x) and salmon (y); the grey bands are changes under 2-fold. Black genes change at "
    f"least 2-fold with both methods in the same direction ({n_both}); blue only with STAR ({n_star}); green only with salmon "
    f"({n_sal}). Panel b: the same counts for the three comparisons. The patch alone gives 5 genes with both methods, all in "
    "its regions. STAR counts a read for a gene only if it maps to one place and overlaps that gene alone (the GTEx-style "
    "standard); salmon, run here on STAR's alignments, keeps reads that fit several places and shares them between the "
    "matching transcripts by likelihood.",
    "Each method has its own way of being wrong. STAR changes whenever a read becomes ambiguous on one reference, because it "
    "throws such reads away; salmon changes whenever the transcript models differ between the two assemblies, which they "
    "do somewhat even in the same RefSeq release. A change produced by the counting method shows up in one method only "
    "(in the pilot, S100A6 changed 10-fold with STAR and not at all with salmon). Requiring both keeps the changes that "
    "come from where the reads go, which is the reference effect we want. It is conservative: some STAR-only genes in "
    "duplications are probably real losses that salmon rescues by sharing reads, so the counts are a lower bound.",
    "Every expression and eQTL result from here on uses this rule. Next: what the agreed changes look like across all genes.")

# 8b ------------------------------------------------------------------------------------------------------------
state['sec'] = '1 · Expression'
# numbers from results/07-presentation/data/expression-compare.tsv (prepare-figure-data.py expression_tables)
ec = pd.read_csv(DD / "expression-compare.tsv", sep="\t")
ep, et = ec[ec.pair == "grch38-ffx-vs-grch38"], ec[ec.pair == "chm13v2.0-vs-grch38-ffx"]
etc = et[et.changed]
sd_all, sd_chg = et.region_class.eq("sd").mean(), etc.region_class.eq("sd").mean()
gcl = pd.read_csv(ROOT / "results/06-compare/gene-classes.tsv", sep="\t", usecols=["gene_ncbi", "t2t_sd98_frac"])
dup = etc.merge(gcl, on="gene_ncbi")
new_n, new_dup = int(dup.region_class.eq("t2t_new_sequence").sum()), int((dup.region_class.eq("t2t_new_sequence") & (dup.t2t_sd98_frac >= 0.5)).sum())
dup_n = int(dup.region_class.eq("sd").sum()) + new_dup   # genes mainly in >= 98% SDs, whatever their label
s = slide(); title(s, "Expression barely changes, except in duplicated genes",
                   "Mean expression in 731 people; a change must be ≥ 2-fold in STAR and in salmon")
px, py, pw, ph = picture(s, fig(14), 0.6, 1.5, 12.1, 4.2)
cols = [(0.0, 0.36, [("Patch alone: ", {"bold": True, "color": FFX}),
                     (f"{int(ep.changed.sum())} genes change, all inside FixItFelix regions. The rest of the genome is untouched.", {})]),
        (0.36, 0.36, [("T2T vs patched: ", {"bold": True, "color": T2T}),
                      (f"{len(etc)} of {len(et):,} genes change; {1 - len(etc) / len(et):.1%} stay within 2-fold. "
                       "Dashed lines: 2-fold.", {})]),
        (0.72, 0.28, [("Where: ", {"bold": True, "color": NAVY}),
                      (f"{sd_chg:.0%} of changed genes are in segmental duplications, against {sd_all:.1%} of all genes.", {})])]
for x0, w0, parts in cols:
    text(s, px + x0 * pw, 5.72, w0 * pw - 0.1, 0.85, parts, size=13)
text(s, 0.6, 6.6, 12.1, 0.4, f"Segmental duplication: ≥ 50% of the gene in copies ≥ 98% identical (T2T SD map). Labels follow "
     f"FixItFelix region > T2T-only sequence > duplication; all {new_n} changed T2T-only genes are duplicated too, so "
     f"{dup_n} of {len(etc)} ({dup_n / len(etc):.0%}) lie in duplications.", size=11, italic=True, color=MUTED)
footer(s)
notes3(s,
    "Gene expression averaged over the 731 people. Each dot is a gene, one reference against the other; on the diagonal "
    "means no change, dashed lines mark 2-fold. Panel a: the patch alone moves only 5 genes, all in its own regions. "
    "Panel b: T2T vs patched GRCh38 moves 97 genes at least 2-fold under both STAR and salmon; 99.3% stay within 2-fold. "
    "Panel c: 31% of the changed genes are in segmental duplications, against 2.4% of all genes, and 7% are in sequence "
    "only T2T has. Genes on the bottom or left edge have almost no reads on one reference: their reads went to another copy. "
    "How the labels work: a segmental duplication is a stretch of DNA, usually over 1 kb, present in two or more "
    "near-identical copies. A gene is labelled duplicated when at least half of it lies in copies at least 98% identical "
    "(the T2T consortium SD map, Vollger et al. 2022); at that identity short reads cannot tell the copies apart. Each "
    "gene gets one label, in the order FixItFelix region, then T2T-only sequence, then duplication, then elsewhere. All 7 "
    "changed T2T-only genes are also fully inside high-identity duplications (copies GRCh38 never assembled), so 37 of "
    "the 97 changed genes (38%) are in duplicated sequence. Elsewhere does not mean unique: it still includes older, "
    "more diverged copies such as many processed pseudogenes, and three genes partly (20-50%) in high-identity duplications.",
    "Two messages. First, T2T does not add noise: almost every gene is measured the same, which is why the overall eQTL "
    "map will be stable. Second, it tells us where to look: the changes are concentrated in duplicated genes and new "
    "sequence, mostly outside the patched regions.",
    "Take-home of the expression section: on the same reads, T2T changes the expression of about 100 genes, confirmed by "
    "both counting methods, almost none of them in the patched regions and mostly in duplicated sequence; the rest of the "
    "transcriptome is measured the same. But an average change says nothing yet about eQTLs, which compare people by "
    "genotype. For that we need the same variant in the same people on both references: the genotypes, next.")

# M_GENO --------------------------------------------------------------------------------------------------------
state['sec'] = '2 · Genotypes'
s = slide(); title(s, "Methods: the same variants on both references", "Genotypes called on each reference, then paired")
method_slide(s, [
    ("Two genotype panels", "1000 Genomes called on T2T (SHAPEIT5-phased, 3,202 people) and on GRCh38 (NYGC phased panel)"),
    ("Keep our people and common variants", "1,089 people with RNA-seq; multi-allelic sites split; allele frequency ≥ 1% in the cohort"),
    ("Pair the variants", "liftOver converts T2T positions to GRCh38; a pair needs the same alleles (or swapped) and frequency within 0.02"),
    ("Ancestry covariates", "5 genotype principal components from pruned common SNPs, shared by all references")],
    [("bcftools ", "1.23.1"), ("plink2 ", "2.0.0a7.8"), ("UCSC liftOver ", "v482, hs1ToHg38, minMatch 0.95")],
    "One list of matched variants; each reference keeps its own genotype file in its own coordinates, filtered to the matched set.")
footer(s)
notes3(s,
    "How the genotypes were prepared. The 1000 Genomes projects called variants separately on T2T and on GRCh38. We kept "
    "the people with RNA-seq and variants with allele frequency of at least 1%. liftOver is used only to find each T2T "
    "variant's GRCh38 position; a pair must also have the same alleles and nearly the same frequency, because it is the "
    "same people. No genotype is ever converted between references.",
    "To compare an eQTL between references we must test the same variant in the same people on both. Pairing guarantees "
    "this, so differences later cannot come from the genotypes.",
    "Next: how many variants pair, and what does not.")

# 9 -------------------------------------------------------------------------------------------------------------
state['sec'] = '2 · Genotypes'
s = slide(); title(s, "A fair comparison of genotypes", "Only variants present on both references enter the main test")
def pvar_counts(name):   # SNP = single-base REF and ALT (multiallelic sites already split)
    v = pd.read_csv(ROOT / f"results/05-eqtl/mage/genotypes/{name}.pvar", sep="\t", comment="#", header=None, usecols=[3, 4])
    snp = (v[3].str.len() == 1) & (v[4].str.len() == 1)
    return len(v) / 1e6, snp.sum() / 1e6, (~snp).sum() / 1e6
gt_all, gg_all, gm = pvar_counts("t2t.all"), pvar_counts("grch38.all"), pvar_counts("t2t.matched")
fig_bullets(s, fig(4), [
    ("Genotypes were called twice, ", "once on each reference, by the 1000 Genomes projects"),
    ("Start (MAGE, frequency ≥ 1%): ", f"T2T {gt_all[0]:.1f} M ({gt_all[1]:.1f} M SNPs, {gt_all[2]:.1f} M indels); "
     f"GRCh38 {gg_all[0]:.1f} M ({gg_all[1]:.1f} M SNPs, {gg_all[2]:.1f} M indels)"),
    ("liftOver only pairs them: ", "same position, same alleles, same frequency; each panel keeps its own coordinates"),
    ("Kept for the eQTL scan: ", f"{gm[0]:.1f} M matched variants ({gm[1]:.1f} M SNPs, {gm[2]:.1f} M indels), on both references (blue)"),
    ("So any change in an eQTL ", "comes from the RNA reads, not from the genotypes")], size=14, after=10)
text(s, 0.5, 6.6, 7.2, 0.4, "MAGE, 731 people, chr1–22. Pairing was done on all 1,089 people with RNA-seq, so variants "
     "below 1% across them were never paired.", size=11, italic=True, color=MUTED)
footer(s)
notes3(s,
    "How the two genotype panels were paired. Variants were called independently on each reference by the 1000 Genomes "
    "projects. liftOver converts positions only; a pair must also have the same two alleles (or swapped), allele "
    "frequency within 0.02, and genotypes that agree person by person. 12.6 million variants pair cleanly. Pairs with "
    "different frequencies (0.13 M) are paralogous positions and are excluded. Some T2T variants have no GRCh38 "
    "position at all: they lie in sequence only T2T has (hypothesis H1).",
    "To say an eQTL changed, we must test exactly the same variant in exactly the same people on both references. If "
    "the genotypes differed, a change could come from the variant calls rather than the reads. With matched variants, "
    "any difference in the eQTL must come from expression, which is the reference effect we want to measure.",
    "Expression and genotypes are ready on every reference. Next: mapping the eQTLs.")

# M_EQTL --------------------------------------------------------------------------------------------------------
state['sec'] = '3 · eQTL mapping'
s = slide(); title(s, "Methods: mapping cis-eQTLs on each reference", "The same recipe on every reference and counting method")
method_slide(s, [
    ("Expression per gene", "expressed genes (TPM ≥ 0.1 and ≥ 6 reads in ≥ 20% of people); TMM-normalised; ranked to a normal scale, so slopes are in SD"),
    ("Covariates (66)", "5 genotype PCs (shared), sex, 60 expression PCs of each reference"),
    ("Test every variant near every gene", "within 1 Mb of the gene start; linear regression of expression on genotype"),
    ("Call eGenes", "10,000 permutations per gene, beta approximation, q < 0.05; all variant–gene results kept for the comparison")],
    [("tensorQTL ", "1.0.10 (map_cis, map_nominal)"), ("plink2 ", "2.0.0a7.8"), ("Python, pandas, edgeR TMM", "")],
    "For each reference and method: eGenes, each gene's lead variant, and the effect of every variant on every nearby gene (chr1–22, matched variants).")
footer(s)
notes3(s,
    "How eQTLs were mapped. For each reference and each counting method, every expressed gene's expression is normalised "
    "and corrected for 66 covariates, then regressed on every matched variant within 1 Mb. Permutations give each gene "
    "an honest p-value given its number of variants, and q < 0.05 across genes defines an eGene. All variant-gene results "
    "are kept, so the comparison can look up any variant's effect on any gene on either reference.",
    "The recipe is identical everywhere, and the genotype PCs are shared, so the only thing that differs between "
    "references is the expression input.",
    "Next: the overall result of the mapping.")

# 10 ------------------------------------------------------------------------------------------------------------
state['sec'] = '3 · eQTL mapping'
s = slide(); title(s, "The overall map is stable", "About 15,000 eGenes on every reference")
fig_bullets(s, fig(5), [
    ("eGene: ", "a gene with at least one significant eQTL (5% FDR, 10,000 permutations, 66 covariates)"),
    ("About 15,000 ", "on every reference and with both counting methods, as in the MAGE paper"),
    ("The totals barely move, ", "so counting eGenes hides which genes change"),
    ("Next: ", "compare the same variant's effect, gene by gene")], fig_w=6.4)
footer(s)
notes3(s,
    "The number of eGenes, genes with at least one significant eQTL (permutations, q < 0.05), on each reference and "
    "with each counting method: about 15,000 everywhere, in line with the MAGE paper on GRCh38.",
    "This is the sanity check for the whole eQTL pipeline, and it says the reference does not change the overall picture. "
    "But totals hide which genes change: PWP2 is an eGene on all three references, yet its effect is four times too "
    "strong on GRCh38. Counting eGenes is the wrong lens.",
    "So what should we compare? The next slide shows why neither expression levels nor eGene lists are enough.")

# 8c ------------------------------------------------------------------------------------------------------------
state['sec'] = '3 · eQTL mapping'
both = int((etc.eqtl_compared & etc.eqtl_changed).sum()); tested = int(etc.eqtl_compared.sum())
eq = et[et.eqtl_changed]; eq_flat = int((~eq.changed).sum())
s = slide(); title(s, "An expression change is not an eQTL change",
                   "An eQTL asks whether expression differs between genotypes, not how high it is")
px, py, pw, ph = picture(s, fig(15), 0.6, 1.5, 12.1, 4.45)
cols = [(0.0, 0.30, [("Same loss for everyone: ", {"bold": True, "color": NAVY}),
                     ("T2T gives every person ~3× fewer reads, whatever their genotype. The genotype difference is "
                      "unchanged, so the eQTL is too (slope 0.74 vs 0.70).", {})]),
        (0.30, 0.30, [("Loss depends on genotype: ", {"bold": True, "color": NAVY}),
                      ("average expression is the same (21 vs 20 CPM), but which genotype gets the reads differs. "
                       "The eQTL reverses (+0.65 vs −0.25).", {})]),
        (0.60, 0.40, [("Genome-wide: ", {"bold": True, "color": NAVY}),
                      (f"of {tested} genes whose expression changes ≥ 2-fold (and that could be tested), only {both} change "
                       f"their eQTL; {eq_flat} of {len(eq)} eQTL changes come with less than a 2-fold expression change.", {})])]
for x0, w0, parts in cols:
    text(s, px + x0 * pw, 6.0, w0 * pw - 0.1, 1.0, parts, size=12)
footer(s)
notes3(s,
    "Why an expression change is not an eQTL change. An eQTL compares people with different genotypes; before mapping, "
    "each gene is ranked and scaled across the 731 people, so only the order of people matters, not the overall level. "
    "Panel a, PKD1P4-NPIPA8: on T2T everyone has about 3 times fewer reads, whatever the genotype; carriers still have "
    "more than non-carriers, so the eQTL is unchanged (slope 0.74 vs 0.70). Panel b, CCZ1: average expression is the "
    "same (21 vs 20 CPM), but CCZ1 shares reads with its copy CCZ1B and the split depends on genotype, so the eQTL "
    "reverses (+0.65 vs -0.25 SD per allele). Panel c: of 59 genes with a 2-fold expression change, only 15 change "
    "their eQTL, and 89 of 104 eQTL changes come with less than a 2-fold expression change.",
    "This justifies the main method: comparing expression or eGene lists would miss most of the real changes and flag "
    "harmless ones. The only fair comparison is the effect of the same variant on the same gene in the same people.",
    "Step 4 does exactly that. First, how to test such a difference correctly when both effects come from the same people.")

# M_CMP ---------------------------------------------------------------------------------------------------------
state['sec'] = '4 · Comparison'
s = slide(); title(s, "Methods: comparing the same eQTL on two references", "Same gene, same variant, same people")
pe_n = len(pd.read_csv(DD / "paired-effects.tsv", sep="\t", usecols=[0]))
z_bonf = norm.isf(0.05 / pe_n / 2)
lead = pd.read_csv(ROOT / "results/05-eqtl/mage/chm13v2.0/star/matched/cis-eqtl.tsv.gz", sep="\t", usecols=["slope", "qval"])
typ = lead.loc[lead.qval < 0.05, "slope"].abs().median()
method_slide(s, [
    ("One variant per gene", "the gene's lead variant and its matched partner on the other reference, alleles counted the same way"),
    ("Its effect on each reference", "slope on T2T and slope on GRCh38, from the same 731 people and the same covariates"),
    ("Paired test of the difference", "difference = slope T2T − slope GRCh38; its error removes the noise both slopes share, "
                                      "because they come from the same people; z = difference ÷ paired error")],
    [("Project scripts ", "compare-effects.py (variant pairs), paired-test.py, define-gene-sets.py"),
     ("ctyper ", "copy numbers (Ma & Chaisson 2025)"), ("numpy, scipy, pandas", "")],
    "For each pair of references: the genes whose eQTL changes, with direction and region; then checked in Geuvadis "
    "and with copy number.", span=3.3)
card(s, 0.45, 4.95, 7.55, 1.95)
dot(s, 0.6, 5.12, 0.48, NAVY, "4", size=15)
text(s, 1.25, 5.06, 6.6, 0.42, "Call a change: all four must hold", size=17, bold=True, color=NAVY)
bullets(s, 1.25, 5.5, 6.65, 1.4, [
    ("Real: ", f"Bonferroni p < 0.05 over about {round(pe_n, -3):,.0f} genes (|z| > {z_bonf:.1f})"),
    ("Big enough: ", f"|difference| ≥ 0.1 SD (a typical eQTL is {typ:.2f} SD)"),
    ("Not one tool: ", "the same with STAR and with salmon, in the same direction"),
    ("Clean: ", "no annotation, MHC or PAR flag")], size=13, after=2)
footer(s)
notes3(s,
    "How the comparison works, in four steps. One: for each gene take its lead variant and the matched copy of that "
    "variant on the other reference, and count the alleles the same way on both, so the two effects are about the same "
    "thing. Two: measure that variant's effect, the slope, on each reference, with the same 731 people and the same "
    "covariates as the eQTL scan. Three: take the difference between the two slopes and divide it by its error. "
    "Because both slopes come from the same people, most of their noise is shared, and the paired error subtracts that "
    "shared part: the error of the difference is se_A squared plus se_B squared minus twice their covariance, "
    "square-rooted. Example: slopes 0.65 and 0.55 give a difference of 0.10; if its paired error is 0.01, z is 10. "
    f"Four: a gene's eQTL has changed only if all four criteria hold. Real: significant after Bonferroni across the "
    f"{pe_n:,} genes compared, which means |z| above {z_bonf:.1f}. Big enough: the slopes differ by at least 0.1 SD of "
    f"expression per allele copy, about a third of a typical eQTL effect ({typ:.2f} SD, the median lead slope); the "
    "paired test can detect changes of 0.05 SD, and 0.1 is where replication in Geuvadis levels off (91%, the same as "
    "at 0.2). Not one tool: the same result with STAR and with salmon, in the same direction. Clean: genes whose "
    "annotation differs between the references, and MHC and PAR genes, are set aside.",
    "Same gene, same variant, same people: the only thing that differs is the reference, so a significant difference "
    "is an effect of the reference. The same steps are run for the three pairs: T2T vs GRCh38, patched vs GRCh38, and "
    "T2T vs patched GRCh38.",
    "Next: the three comparisons, starting with T2T against standard GRCh38.")

# CNT -----------------------------------------------------------------------------------------------------------
state['sec'] = '4 · Comparison'
n02 = int(((tp.delta.abs() >= 0.2) & (tp.delta_salmon.abs() >= 0.2)).sum())
rep_t = tp[tp.geuvadis_same_sign.notna()]
rep_y, rep_n = int(rep_t.geuvadis_nominal.eq(True).sum()), len(rep_t)
s = slide(); title(s, f"{len(tp)} genes change their eQTL beyond the patch", "Main set: real, at least 0.1 SD, STAR and salmon agree")
fig_bullets(s, fig("fig20-changed-counts"), [
    ("T2T vs GRCh38: ", f"{len(tg)} genes change, of about 15,000 compared; {in_ffx['chm13v2.0-vs-grch38']} in patched regions"),
    ("Sanity check, the patch alone: ", f"{len(pg)} change, all inside its regions, as FixItFelix found"),
    ("Beyond the patch: ", f"against patched GRCh38, {len(tp)} still change, none in patched regions"),
    ("Robust: ", f"{n02} remain with a stricter 0.2 SD floor; {rep_y} of {rep_n} replicate in Geuvadis")], fig_w=6.6)
footer(s)
notes3(s,
    f"The headline numbers. A gene counts as changed if its eQTL passes the four criteria from the methods slide. T2T "
    f"vs standard GRCh38: {len(tg)} genes, of which {in_ffx['chm13v2.0-vs-grch38']} lie in FixItFelix's patched regions. "
    f"The sanity check, patched vs standard GRCh38, where the only difference is the patch: {len(pg)} genes change, all "
    f"inside the patched regions, reproducing FixItFelix (CBS, ICOSLG, U2AF1; PWP2 and GATD3 also change but carry an "
    f"annotation flag). T2T vs patched GRCh38: {len(tp)} genes change, none inside the patch. The eGene totals barely "
    f"move (~15,000 on every reference); the change is in which genes and how strongly.",
    f"This is the central result (H0): the known patch explains a small part of what a complete reference changes. It "
    f"does not depend on the threshold: with a 0.2 SD floor {n02} genes remain, none in patched regions, and the "
    f"changes replicate in 358 independent Geuvadis people ({rep_y} of {rep_n} testable genes).",
    "From here on: T2T vs patched GRCh38, what a complete genome adds beyond the known fix. First, which genes these are.")

# TOP10 ---------------------------------------------------------------------------------------------------------
state['sec'] = '4 · Comparison'
t10 = tp.assign(ad=tp.delta.abs()).sort_values("ad", ascending=False).head(10)
btype = {"protein_coding": "Protein-coding", "transcribed_pseudogene": "Pseudogene", "lncRNA": "lncRNA"}
rtype = {"sd": "Seg. duplication", "other": "Elsewhere", "t2t_new_sequence": "T2T-only", "ffx_patched": "FixItFelix"}
yn = lambda v: "Yes" if v is True or v == 1 else ("No" if v is False or v == 0 else "Not testable")
n_sd10 = int((t10.region_class == "sd").sum())
n_rev10 = int((np.sign(t10.slope_a) != np.sign(t10.slope_b)).sum())
s = slide(); title(s, "The changed genes: top 10 by size of change", f"T2T vs patched GRCh38; {len(tp)} genes in all (workbook)")
picture(s, fig("fig23-top10-table"), 0.45, 1.35, 12.4, 4.0)
bullets(s, 0.6, 5.3, 12.1, 1.6, [
    ("Read it like a GWAS table: ", "rsID, position, REF/ALT, the effect allele and its global frequency, z and P"),
    (f"{n_sd10} of 10 ", f"lie in segmental duplications; {n_rev10} of 10 reverse direction between references"),
    ("What to look for next: ", "what kind of genes these are, why their reads move, and whether it touches disease")],
    size=14, after=3)
footer(s)
notes3(s,
    f"The ten genes with the largest change, T2T vs patched GRCh38. Effect is the eQTL slope at the same variant, in SD "
    f"of expression per allele copy; change is T2T minus patched. All ten pass every criterion with STAR and with salmon. "
    f"{n_sd10} of the ten lie in segmental duplications, and {n_rev10} reverse direction: the same allele raises "
    f"expression on one reference and lowers it on the other. Examples: ASAH2 (+0.77 vs −0.53), SMN2 (much stronger "
    f"on T2T), CCZ1. Replicates = same direction with p < 0.05 in the 358 Geuvadis people not in MAGE. The effect allele is the one the "
    f"effects refer to, on the GRCh38 forward strand; its frequency is from 1000 Genomes phase 3 (all 2,504 people), and "
    f"agrees with MAGE except at STAG3L1 (0.56 vs 0.34 in MAGE), a duplicated region where the old low-coverage calls "
    f"are less reliable.",
    "Putting names and numbers up front shows what the rest of the section explains: these are duplicated genes and "
    "pseudogenes, many medically relevant (SMN2, PMS2CL, NBPF10), and their changes are large and reproducible, not "
    "marginal calls.",
    "Next: is the duplication pattern general across all changed genes?")

# TYPES ---------------------------------------------------------------------------------------------------------
state['sec'] = '4 · Comparison'
pe_all = pd.read_csv(DD / "paired-effects.tsv", sep="\t")
pe_all["ch"] = pe_all.tier.isin(["core", "main"]); pe_all["sd"] = pe_all.region_class == "sd"
r_sd = pe_all[pe_all.sd].ch.mean() * 100; r_rest = pe_all[~pe_all.sd].ch.mean() * 100
k_sd, n_sd = int(pe_all[pe_all.sd].ch.sum()), int(pe_all.sd.sum())
sd_all = pe_all.sd.mean() * 100
sd_ch = (tp.region_class == "sd").mean() * 100
sdc = pd.read_csv(DD / "changed-by-sd.tsv", sep="\t").set_index("class")
n_dup, n_all_dup = int(sdc.loc["sd", "changed"]), int(sdc.loc["sd", "all_compared"])
r_dup = n_dup / n_all_dup * 100
r_uni = sdc.loc["unique", "changed"] / sdc.loc["unique", "all_compared"] * 100
s = slide(); title(s, f"{n_dup} of {len(tp)} changed genes lie in segmental duplications",
                   "T2T segmental-duplication map (Vollger et al. 2022): ≥ 90% identical copies covering ≥ half of the gene")
fig_bullets(s, fig("fig21-sd-share"), [
    ("Segmental duplication: ", "a stretch of genome present in two or more near-identical copies"),
    (f"{n_dup / len(tp) * 100:.0f}% of the changes, ", f"but only {n_all_dup / (n_all_dup + int(sdc.loc['unique', 'all_compared'])) * 100:.0f}% "
     f"of the genes compared: {r_dup:.0f}% of duplicated genes change vs {r_uni:.1f}% elsewhere"),
    ("Why: ", "short reads cannot tell such copies apart; GRCh38 misses or collapses many, T2T has them all"),
    ("Unique sequence: ", "often genes partly duplicated or with a copy elsewhere (ASAH2 with ASAH2B, RPSA)"),
    ("Medically relevant families: ", "SMN1/2, PMS2, NBPF, NCF1C")], fig_w=5.6, size=15)
footer(s)
notes3(s,
    f"Where the changed genes are. {sd_ch:.0f}% of the {len(tp)} lie in segmental duplications, stretches present in two "
    f"or more copies at least 98% identical (T2T map, Vollger et al. 2022), against {sd_all:.1f}% of all "
    f"{len(pe_all):,} genes compared. GRCh38 is missing or collapses many of these copies, and short reads cannot tell "
    "near-identical copies apart, so which copy a read is counted on depends on what the reference contains. The other "
    "half include very highly expressed genes with many processed pseudogene copies, such as ribosomal-protein genes.",
    "This points to a mechanism rather than noise: the changes concentrate exactly where read placement is ambiguous. "
    "It also says where eQTL catalogues built on GRCh38 are least reliable.",
    "Next: one gene pair, in the raw RNA coverage.")

# COV -----------------------------------------------------------------------------------------------------------
state['sec'] = '4 · Comparison'
cvd = pd.read_csv(DD / "coverage-example.tsv", sep="\t")
G1, G2 = cvd.gene.drop_duplicates().tolist()[:2]
tot = cvd.groupby(["ref", "gene", "genotype"]).cpm.sum().unstack("genotype")
def shift(ref, gene):
    r = tot.loc[(ref, gene)]
    return r[2] / r[0]
ng = cvd.drop_duplicates("genotype").set_index("genotype").n
v1 = pd.read_csv(DD / "top10-changed.tsv", sep="\t").set_index("symbol").loc[G1]
ea = v1.effect_allele; oa = v1.alt_grch38 if ea == v1.ref_grch38 else v1.ref_grch38
gt0, gt2 = f"{oa}/{oa}", f"{ea}/{ea}"
s = slide(); title(s, "Seeing it in the reads: they move between copies",
                   f"RNA coverage at {G1} and its paralog {G2}, by genotype")
fig_bullets(s, fig("fig25-coverage"), [
    ("Like a genome browser: ", f"mean RNA coverage at {v1.rsid}: {int(ng[0])} people {oa}/{oa} (grey) vs {int(ng[2])} "
     f"{ea}/{ea} (blue; {ea} is the effect allele); exons joined, 5′ to 3′"),
    ("Patched GRCh38 (top): ", f"{ea}/{ea} vs {oa}/{oa}: the {ea} allele lowers {G1} ({shift('grch38-ffx', G1):.2f}×) and its reads go to {G2} "
     f"({shift('grch38-ffx', G2):.2f}×)"),
    ("T2T (bottom): ", f"the {ea} allele raises {G1} ({shift('chm13v2.0', G1):.1f}×); {G2} barely changes "
     f"({shift('chm13v2.0', G2):.2f}×)"),
    ("Same reads, different copy: ", f"in {oa}/{oa} people, reads GRCh38 puts on {G1} sit on {G2} on T2T (grey bump, bottom right)"),
    ("Data: ", "re-aligned RNA of 347 people, uniquely mapped reads (STAR), both references")], fig_w=7.4, size=14)
footer(s)
notes3(s,
    f"A genome-browser view of the mechanism, for the gene with the largest change ({G1}). For the 347 people whose RNA "
    f"was re-aligned for allele-specific expression, we took the coverage of uniquely mapped reads along the exons of "
    f"{G1} and its paralog {G2}, on patched GRCh38 (top) and on T2T (bottom), in reads per million mapped, averaged over "
    f"the people with 0 copies ({int(ng[0])}) and with 2 copies ({int(ng[2])}) of the eQTL allele. The x-axis is the "
    f"exons joined end to end, introns removed, 5′ to 3′. Total coverage, 2 vs 0 copies: patched GRCh38 {G1} "
    f"{shift('grch38-ffx', G1):.2f}×, {G2} {shift('grch38-ffx', G2):.2f}×; T2T {G1} {shift('chm13v2.0', G1):.2f}×, "
    f"{G2} {shift('chm13v2.0', G2):.2f}×.",
    f"This is the picture behind the numbers. On patched GRCh38, people without the allele have reads on {G1}'s 3′ "
    f"exons; on T2T those reads are not on {G1} at all but on {G2}. So the reference decides which copy the reads are "
    f"counted on, and because that depends on the genotype, the eQTL changes sign. It is visible in the raw coverage, "
    f"not an artefact of the test.",
    "Next: the same gene, person by person, and its eQTL on both references.")

# WHY -----------------------------------------------------------------------------------------------------------
state['sec'] = '4 · Comparison'
pdl = pd.read_csv(DD / "paralog-deltas.tsv", sep="\t")
pdl = pdl[pdl.pair == f"{G1} / {G2}"].pivot(index="sample", columns="role", values="delta")
r_cc = np.corrcoef(pdl.gene, pdl.partner)[0, 1]
cz = tp.set_index("symbol").loc[G1]
s = slide(); title(s, "Reads moving reverses the eQTL", f"{G1} and its paralog {G2}, 731 people, T2T vs patched GRCh38")
fig_bullets(s, fig("fig22-reads-move"), [
    ("a · Each dot is a person: ", f"what {G1} loses on T2T, {G2} gains (r = {r_cc:.2f})"),
    ("The reference decides the copy: ", f"the same reads land on {G1} or on {G2}"),
    ("b · Because the move depends on genotype, ", f"the eQTL reverses: {cz.slope_a:+.2f} on T2T, "
     f"{cz.slope_b:+.2f} on patched GRCh38"),
    ("Same pattern: ", "CCZ1/CCZ1B, PMS2/PMS2CL; pseudogenes (RPSA → RPSA2); genes only T2T has")], fig_w=7.0, size=14)
footer(s)
notes3(s,
    f"The consequence for the eQTL. Panel a: for each of the 731 people, the change in {G1} expression between "
    f"references (T2T minus patched GRCh38) against the change in its paralog {G2}. They mirror each other "
    f"(r = {r_cc:.2f}): reads one reference counts on {G1}, the other counts on {G2}. Panel b: {G1} expression by "
    f"genotype on both references. On patched GRCh38 the allele lowers {G1} (slope {cz.slope_b:+.2f}); on T2T it raises "
    f"it ({cz.slope_a:+.2f}). Because the reads that move depend on the person's genotype, the eQTL reverses.",
    "This turns a list of genes into an explanation: the reference decides which copy a read belongs to, and where that "
    "depends on genotype, the eQTL changes. The same pattern appears for CCZ1/CCZ1B (r = -0.96) and PMS2/PMS2CL, for "
    "processed pseudogenes (RPSA to RPSA2) and for genes only T2T contains. 33 of the 57 core genes have one clear partner.",
    "Section 5: does this change which gene a disease signal points to?")

# GWAS1 ---------------------------------------------------------------------------------------------------------
state['sec'] = '5 · Link to GWAS'
cl = pd.read_csv(DD / "coloc.tsv", sep="\t")
cl_m = cl[cl.tier.isin(["core", "main"])]
one = cl_m[cl_m.call.str.endswith("only", na=False)].sort_values("call")
both_n = int((cl_m.call == "both").sum())
s = slide(); title(s, "The disease link of changed genes depends on the reference",
                   "Colocalisation with six GWAS, run on each reference")
picture(s, fig("fig24-gwas-table"), 0.6, 1.45, 12.1, 3.15)
bullets(s, 0.6, 4.85, 12.1, 2.1, [
    ("Colocalisation: ", "across all variants in the region, do the eQTL and the GWAS signal rise and fall together, "
     "i.e. one shared causal variant (PP4) rather than two nearby ones (PP3)? Run on each reference's eQTLs"),
    (f"{len(one)} changed genes ", f"colocalise on one reference only ({both_n} on both); the GWAS data are the same"),
    ("IFITM3-AS1, hypothyroidism: ", "only on patched GRCh38; on T2T its reads go to a gene GRCh38 lacks"),
    ("16p11.2: ", "one ulcerative colitis signal points to SMG1P6 on GRCh38, to NPIPB13 on T2T"),
    ("In short: ", "the reference can change which gene a disease locus is pinned on, so GWAS gene calls in duplicated "
     "regions should be checked on T2T")], size=13, after=3)
footer(s)
notes3(s,
    f"Disease links of the changed genes. For six GWAS (hypothyroidism, lupus, ulcerative colitis, multiple sclerosis, "
    f"type 2 diabetes, atrial fibrillation) we ran colocalisation separately with each reference's eQTLs: the probability "
    f"that the eQTL and the GWAS signal come from the same variant. The GWAS data are the same; only the expression side "
    f"differs. {len(one)} changed genes colocalise on one reference only, {both_n} on both. SuSiE-coloc, which allows more "
    "than one causal variant, confirms IFITM3-AS1, LOC124903630 and SERF1B. Changed genes are not enriched at GWAS loci "
    "(51% vs 46% of matched genes), so the point is not more disease genes but different ones.",
    "This is the practical consequence: the gene a researcher would nominate for a disease locus can change with the "
    "reference. On T2T, IFITM3-AS1's reads go to LOC124906937, a gene GRCh38 does not have, so its GRCh38 link may be "
    "an artefact of missing sequence.",
    "Which reference is right is the open question; next, the limits.")

# 23 ------------------------------------------------------------------------------------------------------------
state['sec'] = 'Summary'
s = slide(); title(s, "Limits and next steps")
for x, head, items, c in [(0.6, "Caveats", [
        "Which reference is right is not yet known (allele-specific expression running)",
        "salmon alone over-calls changes (transcript sets differ between assemblies), so both methods are required",
        "Single-variant coloc is fragile at some loci (PTGES3); 4 loci untestable with European LD",
        "Unmapped 'other' reads: 0.93% on T2T vs 0.09%, not yet explained",
        "24 of 57 core genes have no single read partner yet"], G38),
    (6.8, "Next", [
        "Allele-specific expression at the changed genes and reversals",
        "Read destination for the 60 further main genes",
        "chrX with reference-specific PAR boundaries",
        "GM12878 regulatory elements (H5)"], T2T)]:
    card(s, x, 1.45, 5.9, 5.4)
    dot(s, x + 0.3, 1.7, 0.4, c)
    text(s, x + 0.85, 1.65, 4.8, 0.5, head, size=22, bold=True, color=NAVY, font=HEAD)
    bullets(s, x + 0.3, 2.35, 5.35, 4.4, items, size=15, after=9)
footer(s)
notes3(s,
    "The limits. Which reference is right is not yet known. salmon alone over-calls changes because RefSeq's transcript "
    "sets differ between assemblies, so both quantifiers are required. Single-variant colocalisation is fragile at some "
    "loci (PTGES3). Unmapped 'other' reads are higher on T2T (0.93% vs 0.09%), not yet explained. 24 of 57 core genes "
    "have no single read partner yet.",
    "Saying this upfront makes the main claims stronger: each caveat has a check running or planned, and none of them "
    "affects the central result that most changes lie outside the patched regions.",
    "Next steps: allele-specific expression, read destination for the remaining genes, chrX, and GM12878 regulatory "
    "elements (H5). Then the conclusions.")

# 24 ------------------------------------------------------------------------------------------------------------
state['sec'] = None
s = slide(); title(s, "Conclusions")
concl = [("1", f"T2T vs GRCh38 changes the eQTLs of ~{len(tg)} genes; the overall map (~15,000 eGenes) is stable"),
         ("2", f"The known patch explains {in_ffx['chm13v2.0-vs-grch38']} of them; against patched GRCh38, {len(tp)} genes still change, none in patched regions ({n02} at 0.2 SD)"),
         ("3", "Most changes are reads moving between paralogs, pseudogenes and genes GRCh38 lacks"),
         ("4", "The changes are real: they replicate in Geuvadis (100 of 110) and are not copy number"),
         ("5", "The reference changes which genes colocalise with disease signals (IFITM3-AS1, 16p11.2)")]
for i, (k, t) in enumerate(concl):
    y = 1.35 + i * 0.95
    dot(s, 0.8, y, 0.62, [T2T, FFX, SD, NEW, G38][i], k, size=18)
    text(s, 1.7, y - 0.08, 10.9, 0.8, t, size=19, color=INK, anchor=MSO_ANCHOR.MIDDLE)
text(s, 0.8, 6.15, 11.8, 0.6, "Take-home: check gene families on a complete reference before trusting their eQTLs",
     size=18, bold=True, color=NAVY)
notes3(s,
    f"The five conclusions. (1) T2T vs GRCh38 changes the eQTLs of about {len(tg)} genes while the overall map of ~15,000 "
    f"eGenes is stable. (2) The known patch explains {in_ffx['chm13v2.0-vs-grch38']} of them; against patched GRCh38, "
    f"{len(tp)} still change, none in patched regions. (3) Most changes are reads moving between paralogs, pseudogenes and "
    "genes GRCh38 lacks. (4) The changes are real: they replicate in Geuvadis (100 of "
    "110) and are not copy number. (5) The reference changes which genes colocalise with disease signals.",
    "The take-home for anyone using eQTL catalogues: for gene families and duplicated regions, check the result on a "
    "complete reference before trusting it. Patching known errors is not enough.",
    "Close here; mention that code and results will be released. Backup slides answer questions on numbers, data "
    "sources, tools and the detailed pipeline.")

# Backup 1 ------------------------------------------------------------------------------------------------------
state['sec'] = 'Backup'
s = slide(); title(s, "Backup: key numbers", "MAGE, 731 individuals, chr1-22, matched variants")
eg = pd.read_csv(DD / "egenes.tsv", sep="\t").pivot(index="reference", columns="method", values="egenes")
rows = [["Reference", "eGenes, STAR", "eGenes, salmon"]] + [[n, f"{eg.loc[r, 'star']:,}", f"{eg.loc[r, 'salmon']:,}"]
        for r, n in [("chm13v2.0", "T2T-CHM13"), ("grch38", "GRCh38"), ("grch38-ffx", "Patched GRCh38")]]
table(s, 0.6, 1.75, 5.6, rows, [2.2, 1.7, 1.7], size=13, row_h=0.55, name="eGene table")
tier = lambda p, k: int(sets[p].tier.isin(k).sum())
rows = [["Genes changed", "Core", "Main", "Extent", "In patch"]] + [
    [n, tier(p, ["core"]), tier(p, ["core", "main"]), len(sets[p]), in_ffx[p]] for p, n in PAIRS]
table(s, 6.6, 1.75, 6.1, rows, [2.5, 0.9, 0.9, 0.9, 0.9], size=13, row_h=0.55, name="Change table")
text(s, 0.6, 4.2, 12.1, 0.4, "STAR and salmon agree, no annotation/MHC/PAR flag. Core: independent test; extent: paired test; main: paired, |Δ| ≥ 0.1 SD. In patch: main-set genes in FixItFelix regions.",
     size=12, italic=True, color=MUTED)
for i, (big, lab) in enumerate([("779", "MAGE runs, each aligned to 3 references"), ("6,157", "CPU hours for MAGE expression"),
                                  ("12.6 M", "variants matched between panels"), ("~44", "CPU hours per eQTL scan (reference x method)")]):
    stat(s, 0.6 + i * 3.1, 4.85, 2.8, big, lab, color=NAVY, big_size=34)
footer(s)
notes(s, "Numbers for questions. Full tables are in the workbook (Key results, Gene sets).")

# PWP2 ----------------------------------------------------------------------------------------------------------
state['sec'] = 'Backup'
# numbers from docs/data-journey (STAR counts in each run's ReadsPerGene.out.tab; medians of CPM by genotype)
s = slide(); title(s, "Backup: one gene, two people (PWP2)", "How a reference creates an eQTL · STAR read pairs in MAGE, chr21")
rows = [["Person", "Genotype", "T2T", "GRCh38", "Patched"],
        ["HG00096 (GBR)", "A/A", "613", "11", "628"],
        ["NA19717 (MXL)", "C/C", "828", "313", "851"]]
table(s, 0.6, 1.8, 6.6, rows, [2.1, 1.2, 1.1, 1.1, 1.1], size=16, row_h=0.6, name="PWP2 reads")
text(s, 0.6, 3.75, 6.6, 1.2, [("On GRCh38, HG00096 loses 98% of his reads, NA19717 62%. ", {"bold": True, "color": NAVY}),
     ("A false chr21 copy takes them, and how many it takes depends on the genotype.", {})], size=16)
rows = [["Copies of C", "People", "T2T", "GRCh38", "Patched"],
        ["0", "288", "46.4", "1.3", "48.4"], ["1", "303", "50.0", "10.2", "52.1"], ["2", "140", "50.8", "18.0", "52.8"]]
table(s, 7.6, 1.8, 5.1, rows, [1.3, 0.95, 0.95, 0.95, 0.95], size=15, row_h=0.55, name="PWP2 by genotype")
text(s, 7.6, 4.1, 5.1, 0.4, "Median CPM by genotype at the lead variant", size=12, italic=True, color=MUTED)
for i, (big, lab, c) in enumerate([("0.29", "SD per allele on T2T", T2T), ("1.16", "on GRCh38: four times too strong", G38),
                                    ("0.28", "on patched GRCh38", FFX)]):
    stat(s, 0.6 + i * 4.1, 4.95, 3.8, big, lab, color=c, big_size=40)
footer(s)
notes3(s,
    "PWP2 on chr21, the known artefact, read by read. HG00096 (no C allele) has 613 read pairs on T2T but only 11 on "
    "GRCh38: 98% lost. NA19717 (two C alleles) loses 62%. A false copy of this region on GRCh38 chr21 takes the reads, "
    "and how many it takes depends on the genotype. On T2T and on the patched GRCh38 PWP2 rises gently with each C "
    "allele (0.29 SD per allele); on GRCh38 the slope is 1.16, four times too strong, with q = 1e-210.",
    "This is the whole mechanism in one gene: when the reference makes a person's read loss depend on their genotype, "
    "the eQTL changes. A uniform loss would not matter; a genotype-dependent one creates or distorts an eQTL. And PWP2 "
    "is an eGene on every reference, so simply counting eGenes would miss it.",
    "That is why genotypes come next: to compare an eQTL between references we need the same variant, measured in the "
    "same people, on both.")

# Backup 1b -----------------------------------------------------------------------------------------------------
state['sec'] = 'Backup'
s = slide(); title(s, "Backup: the question and four mechanisms")
card(s, 0.6, 1.45, 12.1, 1.75)
text(s, 0.95, 1.6, 11.4, 1.5, [("H0  ", {"bold": True, "color": NAVY, "size": 26}),
     ("Most reference-driven eQTL changes lie outside the regions FixItFelix patched", {"color": INK, "size": 24})],
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
text(s, 0.6, 6.15, 12.1, 0.7, [("Plan: ", {"bold": True, "color": NAVY}),
     ("compare T2T with GRCh38, check the known patch on its own, then subtract it: T2T vs patched GRCh38.", {})],
     size=16)
footer(s)
notes(s, "The central claim is H0. H1 to H4 are the mechanisms by which a reference can change an eQTL. The talk follows "
         "three comparisons: T2T against GRCh38, the patch alone as a sanity check, and T2T against the patched GRCh38, "
         "which is what is new.")

# PIPE ----------------------------------------------------------------------------------------------------------
state['sec'] = 'Backup'
s = slide(); title(s, "Backup: the detailed pipeline")
picture(s, fig(1), 0.4, 1.15, 7.3, 5.85, align="left")
steps = [("1", "Download and verify", "SRA reads, checksums, read counts"),
         ("2", "Quantify three times", "STAR gene counts and salmon, per reference"),
         ("3", "Genotypes", "Native panels; 12.6 M variants matched"),
         ("4", "Map eQTLs", "tensorQTL, 66 covariates, 10,000 permutations"),
         ("5", "Compare", "Same gene, same variant, paired test")]
for i, (k, h, d) in enumerate(steps):
    y = 1.45 + i * 1.07
    dot(s, 7.95, y, 0.5, NAVY, k, size=15)
    text(s, 8.6, y - 0.07, 4.2, 0.42, h, size=17, bold=True, color=NAVY)
    text(s, 8.6, y + 0.33, 4.2, 0.6, d, size=14, color=MUTED)
footer(s)
notes(s, "Each run is downloaded once and aligned three times, so the only thing that changes is the reference. "
         "STAR counts read pairs placed uniquely on a gene's exons; salmon shares ambiguous reads between transcripts by "
         "likelihood. A change must appear with both. Covariates: 5 genotype PCs shared by all references, sex and 60 "
         "expression PCs. How each tool works is in the backup.")

# Backup 2 ------------------------------------------------------------------------------------------------------
state['sec'] = 'Backup'
s = slide(); title(s, "Backup: where the data come from", "Public data only; every download verified by checksum or read count")
rows = [["Data", "Source / accession", "Publication", "Used for"],
        ["MAGE RNA-seq", "NCBI SRA PRJNA851328", "Taylor 2024, Nature", "Discovery expression"],
        ["Geuvadis RNA-seq", "ENA PRJEB3366", "Lappalainen 2013, Nature", "Replication expression"],
        ["1KGP genotypes, T2T", "Zenodo 7612953 (SHAPEIT5 v1.1)", "Aganezov 2022; Rhie 2023", "T2T-native genotypes"],
        ["1KGP genotypes, GRCh38", "NYGC 2022 phased panel", "Byrska-Bishop 2022, Cell", "GRCh38 genotypes"],
        ["T2T-CHM13v2.0 + RefSeq", "NCBI GCF_009914755.1", "Nurk 2022, Science", "Reference 1"],
        ["GRCh38 no-alt + RefSeq", "NCBI GCA_000001405.15", "Genome Reference Consortium", "Reference 2"],
        ["FixItFelix patched GRCh38", "Zenodo 7532965", "Behera 2023, Genome Biol", "Reference 3, control"],
        ["T2T-only sequence, SDs", "T2T consortium; UCSC chains", "Nurk 2022; Vollger 2022", "Gene region classes"],
        ["ctyper copy numbers", "Zenodo 13683015", "Ma & Chaisson 2025, Nat Genet", "Copy-number covariate"],
        ["GWAS Catalog + 6 GWAS", "EBI GWAS Catalog (GCST)", "Sollis 2023, NAR", "GWAS lookup, coloc"]]
table(s, 0.6, 1.65, 12.1, rows, [3.0, 3.6, 3.2, 2.3], size=12, row_h=0.43, name="Data sources table")
footer(s)
notes(s, "Everything is public. Reads come from NCBI SRA and are checked with vdb-validate and against ENA's read counts. "
         "Gene models are the same NCBI RefSeq release on both assemblies, so genes match by GeneID.")

# Backup 3 ------------------------------------------------------------------------------------------------------
state['sec'] = 'Backup'
s = slide(); title(s, "Backup: how the tools work", "Details and real file examples: docs/data-journey (Quarto)")
tools = [("STAR", "Spliced aligner. A read pair counts for a gene only if it maps to one place and overlaps that gene's exons alone; "
                  "multi-mapped pairs are discarded.", T2T),
         ("salmon", "Shares each fragment between compatible transcripts by likelihood (EM); gives counts and TPM. Keeps multi-mapped reads.", NEW),
         ("Phenotypes", "TPM ≥ 0.1 and ≥ 6 reads in ≥ 20%; TMM-normalised CPM; inverse-normal per gene, so slopes are in SD.", MUTED),
         ("tensorQTL", "Regression per gene × variant (±1 Mb) on 66 covariates; 10,000 permutations + beta fit → gene-level p; q < 0.05 = eGene.", FFX),
         ("Variant pairs", "liftOver positions only; same alleles (or swapped) and MAF within 0.02; sign from genotype correlation.", SD),
         ("Paired test", "z = Δ / √(se_A² + se_B² − 2 cov); cov from the two regressions' residuals (same people).", G38)]
for i, (h, d, c) in enumerate(tools):
    x, y = 0.6 + (i % 2) * 6.15, 1.7 + (i // 2) * 1.75
    card(s, x, y, 5.95, 1.55)
    dot(s, x + 0.25, y + 0.25, 0.32, c)
    text(s, x + 0.7, y + 0.13, 5.0, 0.42, h, size=17, bold=True, color=NAVY)
    text(s, x + 0.7, y + 0.55, 5.1, 0.95, d, size=13, color=INK)
footer(s)
notes(s, "For questions about methods. HG00096's PWP2: 613 STAR pairs, 424 salmon reads, 13.7 TPM on T2T.")

prs.core_properties.title = "A complete reference genome changes the eQTL map beyond known GRCh38 fixes"
prs.core_properties.author = "Yusuf Maulana"
OUT.parent.mkdir(parents=True, exist_ok=True)
# 19 ------------------------------------------------------------------------------------------------------------
state['sec'] = 'Backup'
s = slide(); title(s, "Backup: replication in Geuvadis", "358 individuals not in MAGE, other labs, 2x76 bp reads; the same variants")
gs = pd.read_csv(DD / "gene-sets.tsv", sep="\t")
gm = gs[gs.tier.isin(["core", "main"]) & gs.geuvadis_same_sign.notna()]
h0 = pd.read_csv(DD / "replication-h0.tsv", sep="\t").query("tested")
h1 = pd.read_csv(DD / "replication-h1.tsv", sep="\t")
h1s = h1[(h1.disc_genotype_call == "solid") & h1.disc_signal.isin(["T2T-only variant explains it", "both contribute"])]
rev = h0[h0.effect_class == "Reversed"]
fig_bullets(s, fig(12), [
    ("Independent people: ", "358 Geuvadis individuals not in MAGE, other labs, shorter reads"),
    ("The same test at the same variant pairs ", "as in MAGE"),
    (f"{int(gm.geuvadis_nominal.eq(True).sum())} of {len(gm)} changed genes replicate ",
     f"({int(gm.geuvadis_bonferroni.eq(True).sum())} after Bonferroni)"),
    (f"{int(reversed_again(rev).sum())} of {len(rev)} reversals ", "reverse again"),
    ("H1: ", f"{int(h1s.signal_replicated.sum())} of {len(h1s)} solid T2T-only variant eQTLs replicate")], fig_w=5.9)
footer(s)
notes3(s,
    "Replication in Geuvadis. 93 people are in both cohorts, so we dropped them and used the other 358. For each "
    "main-set gene we took MAGE's exact variant pair and ran the paired test in Geuvadis. 100 of 110 testable genes go "
    "the same way with p < 0.05 under both quantifiers, 86 after Bonferroni; 10 of 12 testable reversals reverse again. "
    "In Geuvadis on its own, T2T vs patched GRCh38 changes 64 genes, again none in FixItFelix regions. For H1, 8 of 15 "
    "well-genotyped T2T-only variant eQTLs replicate, at 5 distinct variants.",
    "Different people, labs and read length give the same answer, so the changes are properties of the references, not "
    "of one dataset. Testing MAGE's exact variant pair means the result does not depend on which lead variant each "
    "cohort happens to pick. H1 contributes examples, not a genome-wide effect.",
    "Next: could something other than the reference, copy number, explain the changes?")

# 20 ------------------------------------------------------------------------------------------------------------
state['sec'] = 'Backup'
s = slide(); title(s, "Backup: copy number does not explain the changes", "ctyper paralog-specific copy numbers as covariates (Ma & Chaisson 2025)")
cn = pd.read_csv(DD / "copy-number.tsv", sep="\t")
cm = cn[cn.tier.isin(["core", "main"])]
fig_bullets(s, fig(18), [
    ("Copy number varies ", "between people in duplicated genes and can drive their expression"),
    ("Added as a covariate, ", "it barely moves the changes: points stay on the diagonal"),
    (f"{int(cm.still_changed.sum())} of {len(cm)} covered genes ", "stay changed; all covered reversals stay"),
    ("Lost (orange): ", ", ".join(cm.loc[~cm.still_changed, "symbol"])),
    ("Caveat: ", f"ctyper covers {len(cm)} of the {len(tp)} genes")], fig_w=5.6)
footer(s)
notes3(s,
    "Copy number as an alternative explanation. Duplicated genes vary in copy number between people, and copy number "
    "can drive their expression (ctyper, Ma and Chaisson 2025). If our lead variants simply tagged copy number, the "
    "references might differ for that reason. Adding gene and family copy numbers as covariates changes almost nothing: "
    "55 of the 62 covered main-set genes keep their change, the median ratio of the change with and without copy number "
    "is 1.00, and all 10 covered reversals stay.",
    "This rules out the most obvious confounder for duplicated genes. The few genes lost (e.g. NUTM2D, whose lead "
    "variant tags copy number) are listed. The caveat: ctyper covers 62 of the 117 genes and counts gene copies, not "
    "every structural variant.",
    "The changes are real and not copy number. Section 5: do they change disease links?")

# 18 ------------------------------------------------------------------------------------------------------------
state['sec'] = 'Backup'
s = slide(); title(s, "Backup: same people, same variant, different answer", "Expression by genotype on two references")
px, py, pw, ph = picture(s, fig(11), 0.6, 1.55, 12.1, 4.5)
cols = [("PWP2: ", "the known artefact; strong only on unpatched GRCh38"), ("CCZ1: ", "the effect reverses direction"),
        ("SMN2: ", "much stronger on T2T")]
for i, (g, d) in enumerate(cols):
    text(s, px + i * pw / 3, 6.15, pw / 3 - 0.1, 0.8, [(g, {"bold": True, "color": NAVY}), (d, {})], size=15)
footer(s)
notes3(s,
    "Expression by genotype for three genes on two references. PWP2: the known artefact, strong only on unpatched GRCh38. "
    "CCZ1: the effect reverses direction. SMN2: much stronger on T2T. 13 core genes reverse: 4 with the variant in "
    "transcribed sequence, where the allele itself steers reads between copies, and 9 with the variant outside exons, "
    "often far away in duplications, probably tagging a structural difference.",
    "The reference can create, remove or reverse an eQTL. Direction is not a detail: colocalisation, Mendelian "
    "randomisation and TWAS all use it, so a reversed eQTL can flip a conclusion about whether a gene raises or "
    "lowers disease risk.",
    "Are these changes real? First, an independent cohort.")

# 22 ------------------------------------------------------------------------------------------------------------
state['sec'] = 'Backup'
s = slide(); title(s, "Backup: which reference is right? Allele-specific expression", "Running now: within-person test of the eQTL direction")
card(s, 0.6, 1.6, 5.9, 5.2)
text(s, 0.9, 1.75, 5.4, 0.45, "The idea", size=18, bold=True, color=NAVY)
bullets(s, 0.9, 2.3, 5.4, 4.4, [
    "In a person heterozygous at the eQTL variant, the haplotype with the raising allele should give more RNA reads",
    "Compared within a person, so covariates and trans effects cancel",
    "Counted on T2T and on patched GRCh38: the reference whose direction matches the allelic imbalance is right"], size=15, after=10)
card(s, 6.8, 1.6, 5.9, 5.2)
text(s, 7.1, 1.75, 5.4, 0.45, "Status", size=18, bold=True, color=NAVY)
bullets(s, 7.1, 2.3, 5.4, 4.4, [
    ("59 of 117 ", "genes testable from short-read SNPs (incl. reversals ELOB, TCAF2, RPS2, PMS2CL)"),
    ("44 SD genes ", "have no exon SNP in the short-read panel, so SNPs are called from 1000 Genomes ONT long reads (260 people)"),
    ("Running: ", "realignment of 347 individuals' RNA in the target genes; long-read calling and phasing (Clair3, WhatsHap)")], size=15, after=10)
footer(s)
notes3(s,
    "Allele-specific expression, running now. In a person heterozygous at the eQTL variant, the haplotype carrying the "
    "raising allele should produce more RNA reads. Because the two alleles are compared inside the same person, "
    "covariates and trans effects cancel. Counting this on T2T and on patched GRCh38 shows which reference's eQTL "
    "direction matches the person's own allelic imbalance. 59 of 117 genes are testable from short-read SNPs; for 44 "
    "duplicated genes SNPs are being called from 1000 Genomes long reads (260 people).",
    "Everything so far shows that the references disagree, not which one is right. This test answers that, and it is "
    "the step that turns 'the reference matters' into 'use T2T for these genes'.",
    "Backup; the limits slide gives the status in one line.")

# 15 ------------------------------------------------------------------------------------------------------------
state['sec'] = 'Backup'
s = slide(); title(s, "Backup: why a paired test", "The standard test treats the two effects as independent")
fig_bullets(s, fig(8), [
    ("For each gene: ", "the same variant's effect on T2T and on patched GRCh38; is the difference real?"),
    ("Same people, same reads: ", "noise is shared and cancels in the difference (residual correlation 0.99)"),
    ("The standard test ", "ignores this and is 9–15× too cautious: inside the grey band nothing would pass"),
    ("The paired test ", "is calibrated (bootstrap 0.95×, null λ = 1.03) and finds the changed genes (blue)")], fig_w=5.9)
footer(s)
notes3(s,
    "The statistics of the comparison. For each gene we take its lead variant, find the same variant on the other "
    "reference, and compare the two slopes. The usual test treats them as independent, but they come from the same "
    "people and mostly the same reads: the median correlation of a gene's residual expression between references is "
    "0.99. Noise that raises a person's expression on one reference raises it on the other too, and cancels in the "
    "difference. The paired test subtracts twice the covariance of the two slopes (Zellner 1962; Lin and Sullivan 2009).",
    "Without it the standard error is 9 to 15 times too large and most real changes are missed. It is not "
    "over-confident: resampling people gives a bootstrap SD 0.95 times the paired standard error, and the null pair "
    "(patched vs unpatched GRCh38 outside the patch) gives lambda 1.03. Example: RPL9's standard error drops from 0.020 "
    "to 0.0075, z from -34 to -92. Changes of at least 0.1 SD replicate in Geuvadis 91% of the time.",
    "Backup, for questions about the statistics of the comparison (methods slide 17).")

# PAIRED2 -------------------------------------------------------------------------------------------------------
state['sec'] = 'Backup'
px = pd.read_csv(DD / "paired-example.tsv", sep="\t")
r0 = px.iloc[0]
r_px = np.corrcoef(px.expr_t2t, px.expr_patched)[0, 1]
sd_ratio = px.expr_t2t.std() / (px.expr_t2t - px.expr_patched).std()
s = slide(); title(s, "Backup: how the paired test is calculated",
                   f"{r0.symbol}, 731 people, STAR, T2T vs patched GRCh38")
fig_bullets(s, fig(19), [
    ("Same people: ", f"each person's expression moves together on both references (r = {r_px:.2f}, panel a)"),
    ("The change is in the difference: ", f"T2T − patched falls with genotype, with {sd_ratio:.1f}× less noise than "
     "expression itself (panel b)"),
    ("Test: ", "z = Δ / √(se²T2T + se²GRCh38 − 2·cov); cov = [Σ eT2T·eGRCh38 / dof] · Σ gT2T·gGRCh38 / (Σ g²T2T · Σ g²GRCh38)"),
    (f"{r0.symbol}: ", f"Δ = {r0.slope_t2t - r0.slope_patched:+.2f} SD; error {r0.se_indep:.3f} if independent, "
     f"{r0.se_paired:.4f} paired; z {r0.z_indep:.1f} → {r0.z_paired:.0f}"),
    ("References: ", "Zellner 1962 (seemingly unrelated regressions); Weesie 1999 (Stata suest); Lin & Sullivan 2009")],
    fig_w=7.0, size=13, after=10)
footer(s)
notes3(s,
    f"How the paired test is calculated, on one gene. Panel a: each dot is one of the 731 people; expression of "
    f"{r0.symbol} after the covariates on patched GRCh38 (x) and T2T (y). The dots hug the diagonal (r = {r_px:.2f}): "
    "whatever makes a person high on one reference makes them high on the other, because the reads are mostly the same. "
    "Panel b: the difference per person, T2T minus patched. It falls with the number of allele copies, which is the eQTL "
    f"change, and its spread is {sd_ratio:.1f} times smaller than that of expression, because the shared noise cancels. "
    "The formula does the same thing without forming the difference by hand: the error of delta is the two standard "
    "errors combined minus twice their covariance; the covariance comes from the two regressions' residuals (e) and "
    "genotypes (g) in the same people, with dof = 731 - 2 - 66 covariates.",
    f"For {r0.symbol} the independent test gives z = {r0.z_indep:.1f}, not significant after Bonferroni; the paired "
    f"test gives z = {r0.z_paired:.0f}. The formula is standard (Zellner 1962; Weesie 1999); Lin & Sullivan (2009) use "
    "the same idea for GWAS with overlapping samples. A bootstrap over people confirms the paired error (0.95×).",
    "Backup, for questions about the statistics of the comparison.")

# COVS ----------------------------------------------------------------------------------------------------------
state['sec'] = 'Backup'
cvs = pd.read_csv(DD / "coverage-example.tsv", sep="\t")
gs1, gs2 = cvs.gene.drop_duplicates().tolist()[:2]
ts = cvs.groupby(["ref", "gene", "genotype"]).cpm.sum().unstack("genotype")
ts = ts.div(ts[0], axis=0)
nn = cvs.drop_duplicates("genotype").set_index("genotype").n
vv = pd.read_csv(DD / "top10-changed.tsv", sep="\t").set_index("symbol").loc[gs1]
e_, o_ = vv.effect_allele, (vv.alt_grch38 if vv.effect_allele == vv.ref_grch38 else vv.ref_grch38)
s = slide(); title(s, f"Backup: {gs1} coverage with all three genotypes",
                   f"Supplementary to slide 21: heterozygotes ({o_}/{e_}) included")
fig_bullets(s, fig("fig25s-coverage-all-genotypes"), [
    ("All 347 re-aligned people: ", f"{int(nn[0])} {o_}/{o_}, {int(nn[1])} {o_}/{e_}, {int(nn[2])} {e_}/{e_}"),
    ("Heterozygotes fall in between ", f"in every panel: the change follows the number of {e_} copies"),
    (f"{gs1} on patched GRCh38: ", f"1 : {ts.loc[('grch38-ffx', gs1), 1]:.2f} : {ts.loc[('grch38-ffx', gs1), 2]:.2f}"),
    (f"{gs1} on T2T: ", f"1 : {ts.loc[('chm13v2.0', gs1), 1]:.2f} : {ts.loc[('chm13v2.0', gs1), 2]:.2f}"),
    ("Why 347, not 731: ", "only the people re-aligned for allele-specific expression have kept BAMs")], fig_w=7.6, size=14)
footer(s)
notes3(s,
    f"The same coverage as slide 21, with the {int(nn[1])} heterozygotes added as a third line. Ratios are total exon "
    f"coverage relative to {o_}/{o_}. On patched GRCh38 {gs1} falls with each copy of {e_} and {gs2} rises; on T2T "
    f"{gs1} rises. Heterozygotes sit between the two homozygous groups, as expected for an additive effect. The "
    f"coverage comes from the 347 MAGE people whose RNA was re-aligned for allele-specific expression (260 with 1000 "
    f"Genomes long-read data plus 87 chosen for heterozygosity at target genes); the full-cohort BAMs were not kept.",
    "Shows the main figure's two-group view is not hiding a non-additive pattern, and gives the sample sizes behind it.",
    "Backup, for questions about the coverage figure.")

prs.save(OUT)
print(OUT, state["n"], "slides")
