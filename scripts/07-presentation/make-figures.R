# Presentation figures 2-19 (ggplot2) from results/07-presentation/data/ (prepare-figure-data.py).
#
#   pixi run Rscript scripts/07-presentation/make-figures.R
#
# Writes results/07-presentation/figures/figNN-*.{png,pdf}: every figure in journal style (<= 183 mm, 7 pt, no titles;
# one idea per figure, the slide carries the explanation), via nat / save_nat. Figure 1 (pipeline) is Graphviz: pipeline.dot.
suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(ggrepel); library(patchwork); library(scales)
})
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "../.."))
dd <- file.path(root, "results/07-presentation/data")
fd <- file.path(root, "results/07-presentation/figures")
dir.create(fd, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) fread(file.path(dd, f), sep = "\t")
save <- function(p, name, w, h) ggsave(file.path(fd, name), p, width = w, height = h, dpi = 300, bg = "white")

theme_set(theme_classic(base_size = 13, base_family = "sans") +
          theme(plot.title = element_text(face = "bold", size = 14), plot.subtitle = element_text(colour = "grey30"),
                strip.background = element_blank(), strip.text = element_text(face = "bold"),
                legend.position = "bottom"))
ref_lab <- c("chm13v2.0" = "T2T-CHM13", "grch38" = "GRCh38", "grch38-ffx" = "Patched GRCh38")
ref_col <- c("T2T-CHM13" = "#0072B2", "GRCh38" = "#D55E00", "Patched GRCh38" = "#E69F00")
cls_lab <- c(other = "Other", sd = "Segmental duplication", t2t_new_sequence = "T2T-only sequence",
             ffx_patched = "FixItFelix region")
cls_col <- c("Other" = "#7F7F7F", "Segmental duplication" = "#CC79A7", "T2T-only sequence" = "#009E73",
             "FixItFelix region" = "#E69F00")
sp_col <- c(AFR = "#E69F00", AMR = "#CC79A7", EAS = "#009E73", EUR = "#0072B2", SAS = "#D55E00")

# Journal style (Nature: <= 183 mm wide, 7 pt text, no titles; the slide carries the explanation). Used by every figure.
# Font: Helvetica metrics (TeX Gyre Heros through cairo/fontconfig). PNG for the deck, PDF for the manuscript.
nat <- theme_classic(base_size = 7, base_family = "TeX Gyre Heros") +
  theme(axis.line = element_line(linewidth = 0.25), axis.ticks = element_line(linewidth = 0.25),
        axis.ticks.length = unit(1.2, "pt"), axis.text = element_text(size = 6.5, colour = "black"),
        axis.title = element_text(size = 7), legend.text = element_text(size = 6.5), legend.title = element_text(size = 7),
        legend.key.size = unit(7, "pt"), legend.background = element_blank(), strip.background = element_blank(),
        strip.text = element_text(size = 7, face = "bold"), plot.tag = element_text(size = 8, face = "bold"),
        plot.margin = margin(2, 4, 2, 2))
save_nat <- function(p, name, w_mm, h_mm) {
  ggsave(file.path(fd, paste0(name, ".png")), p, width = w_mm, height = h_mm, units = "mm", dpi = 300, bg = "white",
         device = png, type = "cairo")
  ggsave(file.path(fd, paste0(name, ".pdf")), p, width = w_mm, height = h_mm, units = "mm", device = cairo_pdf)
}

# Fig 3: mapping QC (journal style) --------------------------------------------------------------------------
qc <- rd("mapping-qc.tsv")
qm <- melt(qc, id.vars = c("run", "sample"), measure.vars = patterns(unique = "^unique_pct_", multi = "^multi_pct_"),
           variable.name = "ref")
qm[, ref := factor(ref_lab[c("chm13v2.0", "grch38", "grch38-ffx")[as.integer(ref)]], levels = ref_lab[c("grch38", "grch38-ffx", "chm13v2.0")])]
map_panel <- function(v, ylab) {
  d <- qm[, .(run, ref, y = get(v))]
  med <- d[, .(m = median(y)), by = ref]
  ggplot(d, aes(ref, y)) +
    geom_point(aes(colour = ref), position = position_jitter(width = 0.17, height = 0, seed = 1), size = 0.35, stroke = 0, alpha = 0.45) +
    geom_boxplot(aes(colour = ref), width = 0.5, outlier.shape = NA, fill = NA, linewidth = 0.35) +
    geom_text(data = med, aes(y = m, label = sprintf("%.1f", m)), nudge_x = 0.36, hjust = 0, size = 6.5 / .pt,
              family = "TeX Gyre Heros") +
    scale_colour_manual(values = ref_col, guide = "none") +
    scale_x_discrete(labels = c("GRCh38" = "GRCh38", "Patched GRCh38" = "Patched\nGRCh38", "T2T-CHM13" = "T2T-CHM13"), expand = expansion(add = c(0.5, 0.85))) +
    labs(x = NULL, y = ylab) + nat
}
p3 <- map_panel("unique", "Uniquely mapped reads (%)") + map_panel("multi", "Multi-mapped reads (%)") +
  plot_annotation(tag_levels = "a")
save_nat(p3, "fig03-mapping-qc", 120, 62)

# Figs 14-15 ----------------------------------------------------------------------------------------------
reg_lab <- c(ffx_patched = "FixItFelix region", sd = "Segmental duplication", t2t_new_sequence = "T2T-only sequence",
             other = "Elsewhere")
reg_col <- c("FixItFelix region" = "#E69F00", "Segmental duplication" = "#CC79A7", "T2T-only sequence" = "#009E73",
             "Elsewhere" = "#000000")
log_x <- scale_x_log10(breaks = 10^(0:4), labels = c("1", "10", "100", "1,000", "10,000"), expand = expansion(mult = 0.02))
log_y <- scale_y_log10(breaks = 10^(0:4), labels = c("1", "10", "100", "1,000", "10,000"), expand = expansion(mult = 0.02))

# Journal-style versions of Figs 2 and 4-13, plus Figs 17-18 (2026-11-02). One idea per figure, no titles; the
# slide carries the explanation. Uses nat / save_nat, ref_col, reg_lab / reg_col (defined above).
lab_txt <- function(...) geom_text(..., size = 6.5 / .pt, family = "TeX Gyre Heros")
repel <- function(...) geom_text_repel(..., size = 6.5 / .pt, family = "TeX Gyre Heros", segment.size = 0.2,
                                       min.segment.length = 0, max.overlaps = 40, seed = 1)
GSD <- file.path(root, "results/06-compare/eqtl-comparison/mage/gene-sets")
pairs3 <- c("chm13v2.0-vs-grch38" = "T2T vs\nGRCh38", "grch38-ffx-vs-grch38" = "Patch vs\nGRCh38",
            "chm13v2.0-vs-grch38-ffx" = "T2T vs\npatched")

# Fig 2: the two cohorts by super-population -------------------------------------------------------------------
co2 <- unique(rd("cohorts.tsv"), by = c("cohort", "sample"))
co2[is.na(superpop) | superpop == "", superpop := "Not given"]
co2[, cohort := factor(fifelse(cohort == "mage", "MAGE\n(discovery)", "Geuvadis\n(replication)"),
                       levels = c("Geuvadis\n(replication)", "MAGE\n(discovery)"))]
c2 <- co2[, .N, by = .(cohort, superpop)]
p2 <- ggplot(c2, aes(N, cohort, fill = superpop)) + geom_col(width = 0.6) +
  geom_text(data = c2[N >= 40], aes(label = N), position = position_stack(vjust = 0.5), colour = "white",
            size = 6.5 / .pt, family = "TeX Gyre Heros") +
  lab_txt(data = c2[, .(N = sum(N)), by = cohort], aes(N, cohort, label = paste0("n = ", N)), inherit.aes = FALSE, hjust = -0.15) +
  scale_fill_manual(values = c(sp_col, "Not given" = "grey70"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Individuals", y = NULL) + nat + theme(legend.position = "bottom")
save_nat(p2, "fig02-cohorts", 100, 45)

# Fig 4: variant pairing, as seen by the MAGE eQTL scan (chr1-22, MAF >= 1% within MAGE's 731 people) ---------------
vm <- rd("variant-matching-mage.tsv")[category != "all", .(n = sum(n)), by = .(category)]
vm_lev <- c(matched = "Matched: same variant on both", matched_maf_discordant = "Paired, frequency differs",
            lifted_not_called_grch38 = "T2T variant, not called on GRCh38", no_liftover = "In T2T-only sequence",
            not_paired_rare_in_1089 = "Not paired (< 1% across all 1,089)",
            grch38_no_matched_partner = "GRCh38 variant, no matched T2T partner")
vm[, label := factor(vm_lev[category], levels = rev(vm_lev))]
vm[, used := category == "matched"]
p4 <- ggplot(vm, aes(n / 1e6, label, fill = used)) + geom_col(width = 0.65) +
  lab_txt(aes(label = sprintf("%.2f M", n / 1e6)), hjust = -0.12) +
  scale_fill_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "grey70"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(x = "MAGE variants, allele frequency \u2265 1%, chr1\u201322 (millions)", y = NULL) + nat
save_nat(p4, "fig04-variant-matching", 110, 50)

# Fig 5: eGenes per reference -------------------------------------------------------------------------------------
eg <- rd("egenes.tsv")
eg[, ref := factor(ref_lab[reference], levels = ref_lab)]
eg[, method := factor(fifelse(method == "star", "STAR", "salmon"), levels = c("STAR", "salmon"))]
p5 <- ggplot(eg, aes(method, egenes, fill = ref)) + geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_text(aes(label = comma(egenes)), position = position_dodge(0.75), vjust = -0.4, size = 6 / .pt, family = "TeX Gyre Heros") +
  scale_fill_manual(values = ref_col, name = NULL) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "eGenes (FDR 5%)") + nat + theme(legend.position = "bottom")
save_nat(p5, "fig05-egenes", 80, 60)

# Fig 6: positive control ------------------------------------------------------------------------------------------
pcq <- rd("positive-control.tsv")
pcq[, ref := factor(ref_lab[reference], levels = ref_lab[c("grch38", "grch38-ffx", "chm13v2.0")])]
p6 <- ggplot(pcq, aes(ref, -log10(qval), group = symbol)) +
  geom_hline(yintercept = -log10(0.05), linetype = "22", linewidth = 0.25, colour = "grey60") +
  geom_line(linewidth = 0.35, colour = "grey50") + geom_point(aes(colour = ref), size = 1.4) +
  repel(data = pcq[reference == "grch38"], aes(label = symbol), nudge_x = -0.35, direction = "y", hjust = 1) +
  scale_colour_manual(values = ref_col, guide = "none") +
  scale_y_log10() + scale_x_discrete(expand = expansion(add = c(0.8, 0.3))) +
  labs(x = NULL, y = expression(paste(-log[10], " q (eGene strength)"))) + nat
save_nat(p6, "fig06-positive-control", 80, 65)

# Fig 7: same variant, two references ------------------------------------------------------------------------------
ef <- rd("effects.tsv")
ef[, grp := factor(fcase(tier %in% c("core", "main"), "Changed (main set)", tier == "extent", "Smaller change",
                         default = "No change"), levels = c("Changed (main set)", "Smaller change", "No change"))]
n7 <- ef[, .N, by = grp]
lab7 <- setNames(sprintf("%s (%s)", levels(ef$grp), comma(n7[match(levels(ef$grp), grp), N])), levels(ef$grp))
ef[, absd := abs(slope_a - slope_b)]
show7 <- ef[symbol %in% c("CCZ1", "ASAH2", "SMN2", "RPL9", "NBPF10", "PMS2", "RPSA")]
p7 <- ggplot(ef[order(-as.integer(grp))], aes(slope_b, slope_a, colour = grp)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.25, colour = "grey60") +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey85") + geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey85") +
  geom_point(aes(size = grp), stroke = 0) +
  repel(data = show7, aes(label = symbol), colour = "black", fontface = "italic") +
  scale_size_manual(values = c(1.1, 0.6, 0.35), guide = "none") +
  scale_colour_manual(values = c("Changed (main set)" = "#0072B2", "Smaller change" = "#56B4E9", "No change" = "grey80"),
                      labels = lab7, name = NULL) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 1.6))) +
  coord_equal() +
  labs(x = "Effect on patched GRCh38 (SD per allele)", y = "Effect on T2T-CHM13 (SD per allele)") +
  nat + theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1))
save_nat(p7, "fig07-effect-comparison", 85, 85)

# Fig 8: paired vs independent test --------------------------------------------------------------------------------
pe <- rd("paired-effects.tsv")[artefact_flag == FALSE]
zb <- qnorm(1 - 0.025 / nrow(pe))
plog <- pseudo_log_trans(sigma = 2, base = 10)
zbr <- c(-100, -10, -3, 0, 3, 10, 100)
pe[, grp := factor(fifelse(tier %in% c("core", "main"), "Changed (main set)", "Other genes"), levels = c("Changed (main set)", "Other genes"))]
p8 <- ggplot(pe[order(grp == "Changed (main set)")], aes(z_diff, z_paired, colour = grp)) +
  annotate("rect", xmin = -zb, xmax = zb, ymin = -Inf, ymax = Inf, fill = "grey95") +
  geom_hline(yintercept = c(-zb, zb), linetype = "22", linewidth = 0.25, colour = "grey60") +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.25, colour = "grey70") +
  geom_point(aes(size = grp), stroke = 0) + scale_size_manual(values = c(1, 0.35), guide = "none") +
  scale_colour_manual(values = c("Changed (main set)" = "#0072B2", "Other genes" = "grey75"), name = NULL) +
  scale_x_continuous(trans = plog, breaks = zbr) + scale_y_continuous(trans = plog, breaks = zbr) +
  guides(colour = guide_legend(override.aes = list(size = 1.6))) +
  labs(x = "z, standard test (effects treated as independent)", y = "z, paired test (same people)") +
  nat + theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1))
save_nat(p8, "fig08-paired-test", 85, 80)

# Fig 9: what the 117 changed genes are -----------------------------------------------------------------------------
gs <- rd("gene-sets.tsv")[tier %in% c("core", "main")]
gs[, pattern := factor(fcase(sign(slope_a) != sign(slope_b), "Opposite\ndirection", abs(slope_a) > abs(slope_b), "Stronger\non T2T",
                             default = "Stronger on\npatched GRCh38"),
                       levels = c("Stronger\non T2T", "Stronger on\npatched GRCh38", "Opposite\ndirection"))]
gs[, region := factor(reg_lab[region_class], levels = reg_lab)]
c9 <- gs[, .N, by = .(pattern, region)]
p9 <- ggplot(c9, aes(pattern, N, fill = region)) + geom_col(width = 0.65) +
  geom_text(data = c9[N >= 4], aes(label = N), position = position_stack(vjust = 0.5), colour = "white",
            size = 6.5 / .pt, family = "TeX Gyre Heros") +
  lab_txt(data = c9[, .(N = sum(N)), by = pattern], aes(pattern, N, label = N), inherit.aes = FALSE, vjust = -0.4) +
  scale_fill_manual(values = reg_col, drop = TRUE, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Changed genes (main set)") + nat + theme(legend.position = "right")
save_nat(p9, "fig09-robust-genes", 100, 60)

# Fig 10: reads move between copies ----------------------------------------------------------------------------------
pd <- rd("paralog-deltas.tsv")
pw <- dcast(pd, pair + sample ~ role, value.var = "delta")
r10 <- pw[, .(r = cor(gene, partner)), by = pair]
pw <- merge(pw, r10, by = "pair")
pw[, pair := factor(sprintf("%s  (r = %.2f)", pair, r))]
p10 <- ggplot(pw, aes(gene, partner)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey85") + geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey85") +
  geom_abline(slope = -1, intercept = 0, linetype = "22", linewidth = 0.25, colour = "grey50") +
  geom_point(size = 0.4, stroke = 0, alpha = 0.6, colour = "#CC79A7") +
  facet_wrap(~pair, scales = "free", ncol = 2) +
  labs(x = "Change in the gene, T2T − patched GRCh38 (CPM per person)", y = "Change in its partner (CPM)") +
  nat + theme(strip.text = element_text(size = 6.5, face = "bold.italic"))
save_nat(p10, "fig10-paralog-swaps", 100, 85)

# Fig 11: the same variant, a different eQTL --------------------------------------------------------------------------
ex11 <- rd("eqtl-examples.tsv")
ex11[, ref := factor(ref_lab[reference], levels = ref_lab)]
ex11[, g := factor(round(genotype), levels = 0:2)]
ex11[, symbol := factor(symbol, levels = c("PWP2", "CCZ1", "SMN2"))]
p11 <- ggplot(ex11[!is.na(g)], aes(g, expression, colour = ref)) +
  geom_boxplot(position = position_dodge(0.75), width = 0.6, outlier.shape = NA, fill = NA, linewidth = 0.3) +
  facet_wrap(~symbol, nrow = 1) + scale_colour_manual(values = ref_col, name = NULL) +
  labs(x = "Copies of the eQTL allele", y = "Normalised expression (SD)") +
  nat + theme(legend.position = "bottom", strip.text = element_text(face = "bold.italic"))
save_nat(p11, "fig11-eqtl-examples", 150, 62)

# Fig 12: replication in Geuvadis -------------------------------------------------------------------------------------
rp12 <- rd("replication-paired.tsv")[tier %in% c("core", "main") & !is.na(z_paired)]
rp12[, ok := factor(fifelse(nominal, "Replicates (p < 0.05, same sign)", "Does not replicate"),
                    levels = c("Replicates (p < 0.05, same sign)", "Does not replicate"))]
n12 <- rp12[, .N, by = ok]
lab12 <- setNames(sprintf("%s: %d", levels(rp12$ok), n12[match(levels(rp12$ok), ok), N]), levels(rp12$ok))
p12 <- ggplot(rp12, aes(mage_z_paired, z_paired, colour = ok)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey85") + geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey85") +
  geom_abline(slope = 1, intercept = 0, linetype = "22", linewidth = 0.25, colour = "grey60") +
  geom_point(size = 1, stroke = 0) +
  repel(data = rp12[ok == "Does not replicate"], aes(label = symbol), colour = "grey40", fontface = "italic") +
  scale_colour_manual(values = c("Replicates (p < 0.05, same sign)" = "#0072B2", "Does not replicate" = "grey60"),
                      labels = lab12, name = NULL) +
  scale_x_continuous(trans = plog, breaks = zbr) + scale_y_continuous(trans = plog, breaks = zbr) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 1.6))) +
  labs(x = "Paired z, MAGE (discovery)", y = "Paired z, Geuvadis (358 people not in MAGE)") +
  nat + theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1))
save_nat(p12, "fig12-replication", 85, 80)

# Fig 13: colocalisation with GWAS on each reference ------------------------------------------------------------------
co13 <- rd("coloc.tsv")
setnames(co13, c("PP4_chm13v2.0", "PP4_grch38-ffx"), c("pp4_t2t", "pp4_ffx"))
co13[, cl := factor(fcase(call == "both", "Both references", call == "chm13v2.0 only", "T2T only",
                          call == "grch38-ffx only", "Patched GRCh38 only", default = "Suggestive"),
                    levels = c("T2T only", "Patched GRCh38 only", "Both references", "Suggestive"))]
co13[, lab := paste0(symbol, ", ", fcase(trait == "hypothyroidism", "hypothyroidism", trait == "atrial fibrillation", "AF",
                                          trait == "ulcerative colitis", "UC", trait == "type 2 diabetes", "T2D",
                                          trait == "systemic lupus erythematosus", "SLE", default = trait))]
lo13 <- co13[cl == "Patched GRCh38 only"][order(-pp4_ffx, -pp4_t2t)]
lo13[, ly := 0.42 - 0.075 * (seq_len(.N) - 1)]
p13 <- ggplot(co13, aes(pp4_ffx, pp4_t2t, colour = cl)) +
  annotate("rect", xmin = -0.02, xmax = 0.5, ymin = 0.8, ymax = 1.02, fill = "#0072B2", alpha = 0.08) +
  annotate("rect", xmin = 0.8, xmax = 1.02, ymin = -0.02, ymax = 0.5, fill = "#E69F00", alpha = 0.12) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.25, colour = "grey80") + geom_point(size = 1.3, stroke = 0) +
  repel(data = co13[cl == "T2T only"], aes(label = lab), show.legend = FALSE) +
  geom_segment(data = lo13, aes(x = 0.6, y = ly, xend = pp4_ffx, yend = pp4_t2t), linewidth = 0.2, colour = "grey60",
               inherit.aes = FALSE) +
  geom_text(data = lo13, aes(x = 0.59, y = ly, label = lab), hjust = 1, colour = "#E69F00", size = 6.5 / .pt,
            family = "TeX Gyre Heros", inherit.aes = FALSE) +
  scale_colour_manual(values = c("T2T only" = "#0072B2", "Patched GRCh38 only" = "#E69F00", "Both references" = "grey30",
                                 "Suggestive" = "grey70"), name = NULL) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 1.6))) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Colocalisation probability on patched GRCh38", y = "Colocalisation probability on T2T") +
  nat + theme(legend.position = "right")
save_nat(p13, "fig13-gwas", 120, 80)

# Fig 17: changed genes by region, the three comparisons (main sets) ----------------------------------------------------
cs <- rbindlist(lapply(names(pairs3), function(p) {
  d <- fread(file.path(GSD, paste0("changed-", p, ".tsv")))[tier %in% c("core", "main")]
  d[, .N, by = region_class][, pair := pairs3[[p]]]
}))
cs[, pair := factor(pair, levels = pairs3)]
cs[, region := factor(reg_lab[region_class], levels = reg_lab)]
bar17 <- function(d) ggplot(d, aes(pair, N, fill = region)) + geom_col(width = 0.6) +
  geom_text(data = d[N >= 5], aes(label = N), position = position_stack(vjust = 0.5), colour = "white",
            size = 6.5 / .pt, family = "TeX Gyre Heros") +
  lab_txt(data = d[, .(N = sum(N)), by = pair], aes(pair, N, label = N), inherit.aes = FALSE, vjust = -0.4) +
  scale_fill_manual(values = reg_col, drop = TRUE, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Changed genes (main set)") + nat + theme(legend.position = "right")
save_nat(bar17(cs), "fig17-three-comparisons", 100, 65)
save_nat(bar17(cs[pair == pairs3[["chm13v2.0-vs-grch38"]]]), "fig17a-t2t-vs-grch38", 75, 65)

# Fig 18: copy number as a covariate --------------------------------------------------------------------------------------
cn18 <- rd("copy-number.tsv")[tier %in% c("core", "main")]
cn18[, kept := factor(fifelse(still_changed, "Still changed", "Lost with copy number"), levels = c("Still changed", "Lost with copy number"))]
n18 <- cn18[, .N, by = kept]
lab18 <- setNames(sprintf("%s (%d)", levels(cn18$kept), n18[match(levels(cn18$kept), kept), N]), levels(cn18$kept))
p18 <- ggplot(cn18, aes(delta, delta_cn, colour = kept)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey85") + geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey85") +
  geom_abline(slope = 1, intercept = 0, linetype = "22", linewidth = 0.25, colour = "grey60") +
  geom_point(size = 1.1, stroke = 0) +
  repel(data = cn18[kept == "Lost with copy number"], aes(label = symbol), fontface = "italic", show.legend = FALSE) +
  scale_colour_manual(values = c("Still changed" = "#0072B2", "Lost with copy number" = "#D55E00"), labels = lab18, name = NULL) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 1.6))) + coord_equal() +
  labs(x = "Change in effect, without copy number (SD)", y = "Change in effect, with copy number (SD)") +
  nat + theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1))
save_nat(p18, "fig18-copy-number", 80, 80)


ec <- rd("expression-compare.tsv")
ec[, region := factor(reg_lab[region_class], levels = reg_lab)]
scatter <- function(pr, xl, yl) {
  d <- ec[pair == pr]; d[, `:=`(x = pmax(cpm_b_star, 0.5), y = pmax(cpm_a_star, 0.5))]
  ch <- d[changed == TRUE]
  ggplot(d, aes(x, y)) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.25, colour = "grey60") +
    geom_abline(slope = 1, intercept = c(-1, 1) * log10(2), linewidth = 0.25, colour = "grey60", linetype = "22") +
    geom_point(colour = "grey80", size = 0.25, stroke = 0) +
    geom_point(data = ch, aes(fill = region), shape = 21, size = 1.0, stroke = 0) +
    annotate("text", x = 0.6, y = 2.2e4, hjust = 0, vjust = 1, size = 6.5 / .pt, family = "TeX Gyre Heros",
             label = sprintf("%s genes\n%d change ≥ 2-fold", format(nrow(d), big.mark = ","), nrow(ch))) +
    scale_fill_manual(values = reg_col, drop = FALSE, guide = "none") + log_x + log_y +
    coord_equal(xlim = c(0.5, 2.5e4), ylim = c(0.5, 2.5e4)) +
    labs(x = xl, y = yl) + nat
}
p14a <- scatter("grch38-ffx-vs-grch38", "GRCh38 (mean CPM)", "Patched GRCh38 (mean CPM)")
p14b <- scatter("chm13v2.0-vs-grch38-ffx", "Patched GRCh38 (mean CPM)", "T2T-CHM13 (mean CPM)")
sh <- rbind(ec[pair == "chm13v2.0-vs-grch38-ffx", .(n = .N), by = region][, grp := "All\nexpressed"],
            ec[pair == "chm13v2.0-vs-grch38-ffx" & changed == TRUE, .(n = .N), by = region][, grp := "T2T vs\npatched"],
            ec[pair == "grch38-ffx-vs-grch38" & changed == TRUE, .(n = .N), by = region][, grp := "Patch vs\nGRCh38"])
sh[, `:=`(f = n / sum(n), N = sum(n)), by = grp]
sh[, grp := factor(grp, levels = c("All\nexpressed", "Patch vs\nGRCh38", "T2T vs\npatched"))]
p14c <- ggplot(sh, aes(grp, f, fill = region)) + geom_col(width = 0.65, colour = NA) +
  geom_text(data = unique(sh[, .(grp, N)]), aes(grp, 1.02, label = format(N, big.mark = ",", trim = TRUE)),
            inherit.aes = FALSE, vjust = 0, size = 6.5 / .pt, family = "TeX Gyre Heros") +
  geom_text(data = sh[f >= 0.07], aes(label = percent(f, accuracy = 1)), position = position_stack(vjust = 0.5),
            size = 6.5 / .pt, family = "TeX Gyre Heros", colour = "white") +
  scale_fill_manual(values = reg_col, drop = FALSE, name = NULL) +
  scale_y_continuous(labels = percent, breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Share of genes") + nat
save_nat((p14a + p14b + p14c + plot_layout(widths = c(1, 1, 0.75), guides = "collect") +
          plot_annotation(tag_levels = "a")) & theme(legend.position = "bottom"),
         "fig14-expression", 183, 78)

# Fig 15: an expression change is not an eQTL change ----------------------------------------------------------------
ex <- rd("expression-eqtl-examples.tsv")
ex[, ref := factor(ref_lab[reference], levels = ref_lab[c("chm13v2.0", "grch38-ffx")])]
ex_panel <- function(sym, ylim_top) {
  d <- ex[symbol == sym]
  d[, g := factor(genotype, levels = 0:2)]
  sl <- unique(d[, .(ref, slope)])
  sl[, txt := sprintf("slope %+.2f", slope)]
  ggplot(d, aes(g, cpm, colour = ref)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.75, seed = 1), size = 0.3, stroke = 0, alpha = 0.5) +
    geom_boxplot(position = position_dodge(0.75), width = 0.6, outlier.shape = NA, fill = NA, linewidth = 0.3) +
    geom_text(data = sl, aes(x = 0.5, y = ylim_top * c(1, 0.9)[as.integer(ref)], label = txt), hjust = 0, vjust = 1,
              size = 6.5 / .pt, family = "TeX Gyre Heros", show.legend = FALSE) +
    scale_colour_manual(values = ref_col[levels(ex$ref)], name = NULL) +
    scale_x_discrete(drop = TRUE) + coord_cartesian(ylim = c(0, ylim_top), expand = FALSE, xlim = c(0.4, NA)) +
    annotate("text", x = Inf, y = ylim_top, label = sym, hjust = 1, vjust = 1, size = 7 / .pt, fontface = "italic",
             family = "TeX Gyre Heros") +
    labs(x = "Copies of the eQTL allele", y = "Expression (CPM)") + nat
}
p15a <- ex_panel("PKD1P4-NPIPA8", 50)
p15b <- ex_panel("CCZ1", 60)
pe <- rd("paired-effects.tsv")[, .(gene_ncbi = as.character(gene_ncbi), delta = slope_a - slope_b)]
sc <- merge(ec[pair == "chm13v2.0-vs-grch38-ffx" & eqtl_compared == TRUE][, gene_ncbi := as.character(gene_ncbi)], pe, by = "gene_ncbi")
sc[, cls := factor(fcase(changed & eqtl_changed, "Both", changed, "Expression only", eqtl_changed, "eQTL only", default = "Neither"),
                   levels = c("Expression only", "eQTL only", "Both", "Neither"))]
cnt <- sc[, .N, by = cls]
cls_col2 <- c("Expression only" = "#56B4E9", "eQTL only" = "#D55E00", "Both" = "#000000", "Neither" = "grey80")
lab2 <- setNames(sprintf("%s (%d)", levels(sc$cls), cnt[match(levels(sc$cls), cls), N]), levels(sc$cls))
p15c <- ggplot(sc[order(cls != "Neither")], aes(log2fc_star, delta, colour = cls)) +
  geom_vline(xintercept = c(-1, 1), linewidth = 0.25, linetype = "22", colour = "grey60") +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey60") +
  geom_point(size = 0.5, stroke = 0) +
  geom_point(data = sc[symbol %in% c("PKD1P4-NPIPA8", "CCZ1")], shape = 21, size = 1.8, stroke = 0.4, colour = "black", fill = NA) +
  geom_text(data = sc[symbol %in% c("PKD1P4-NPIPA8", "CCZ1")][, `:=`(hj = fifelse(symbol == "CCZ1", 1.15, 0.3), vj = fifelse(symbol == "CCZ1", -0.5, -1.2))],
            aes(label = symbol, hjust = hj, vjust = vj), colour = "black", size = 6.5 / .pt, family = "TeX Gyre Heros", fontface = "italic") +
  scale_colour_manual(values = cls_col2, labels = lab2, name = "Changed:") +
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 1.5))) +
  coord_cartesian(xlim = c(-4, 4), ylim = c(-1.1, 1.1)) +
  labs(x = expression(paste("Expression change, ", log[2], " T2T / patched GRCh38")),
       y = "eQTL change (slope difference, SD)") + nat
save_nat((p15a + p15b + plot_layout(guides = "collect") & theme(legend.position = "bottom")) - p15c +
          plot_layout(widths = c(2, 1.2)) + plot_annotation(tag_levels = "a") & theme(legend.position = "bottom"),
         "fig15-expression-vs-eqtl", 183, 80)


# Fig 16: STAR vs salmon, the two ways of counting (journal style) -------------------------------------------------
# Each gene's expression change between two references under STAR (unique reads only) and salmon (ambiguous reads
# shared by likelihood). A change counts only when both see >= 2-fold in the same direction.
ec[, call := factor(fcase(changed, "Both methods",
                          abs(log2fc_star) >= 1, "STAR only",
                          abs(log2fc_salmon) >= 1, "salmon only",
                          default = "Neither"), levels = c("Both methods", "STAR only", "salmon only", "Neither"))]
call_col <- c("Both methods" = "#000000", "STAR only" = "#0072B2", "salmon only" = "#009E73", "Neither" = "grey80")
d16 <- ec[pair == "chm13v2.0-vs-grch38-ffx"]
n16 <- d16[, .N, by = call]
lab16 <- setNames(sprintf("%s (%s)", levels(d16$call), format(n16[match(levels(d16$call), call), N], big.mark = ",", trim = TRUE)),
                  levels(d16$call))
clamp <- function(v) pmin(pmax(v, -6), 6)
p16a <- ggplot(d16[order(call != "Neither")], aes(clamp(log2fc_star), clamp(log2fc_salmon), colour = call)) +
  annotate("rect", xmin = -1, xmax = 1, ymin = -Inf, ymax = Inf, fill = "grey95") +
  annotate("rect", ymin = -1, ymax = 1, xmin = -Inf, xmax = Inf, fill = "grey95") +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.25, colour = "grey60") +
  geom_point(size = 0.5, stroke = 0) +
  scale_colour_manual(values = call_col, labels = lab16, name = NULL) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 1.5))) +
  coord_equal(xlim = c(-6, 6), ylim = c(-6, 6)) +
  scale_x_continuous(breaks = seq(-6, 6, 2)) + scale_y_continuous(breaks = seq(-6, 6, 2)) +
  labs(x = expression(paste("STAR: ", log[2], " T2T / patched GRCh38")),
       y = expression(paste("salmon: ", log[2], " T2T / patched GRCh38"))) +
  nat + theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1))
b16 <- ec[call != "Neither", .N, by = .(pair, call)]
b16[, pair := factor(c("chm13v2.0-vs-grch38" = "T2T vs\nGRCh38", "grch38-ffx-vs-grch38" = "Patch vs\nGRCh38",
                       "chm13v2.0-vs-grch38-ffx" = "T2T vs\npatched")[pair],
                     levels = c("T2T vs\nGRCh38", "Patch vs\nGRCh38", "T2T vs\npatched"))]
p16b <- ggplot(b16, aes(pair, N, fill = call)) + geom_col(width = 0.65) +
  geom_text(data = b16[N >= 15], aes(label = N), position = position_stack(vjust = 0.5), size = 6.5 / .pt,
            family = "TeX Gyre Heros", colour = "white") +
  geom_text(data = b16[, .(N = sum(N), both = sum(N[call == "Both methods"])), by = pair],
            aes(pair, N, label = sprintf("%d, %d both", N, both)), inherit.aes = FALSE, vjust = -0.4, size = 6.5 / .pt,
            family = "TeX Gyre Heros") +
  scale_fill_manual(values = call_col, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Genes changed ≥ 2-fold") + nat
save_nat(p16a + p16b + plot_layout(widths = c(1.35, 1)) + plot_annotation(tag_levels = "a"), "fig16-star-salmon", 140, 72)

cat("figures written to", fd, "\n")

# Fig 19: how the paired test works, one gene (IFITM2, STAR, T2T vs patched GRCh38) ---------------------------------
px <- rd("paired-example.tsv")
px[, g := factor(round(genotype), levels = 0:2)]
gcol <- c("0" = "#BBBBBB", "1" = "#56B4E9", "2" = "#0072B2")
r19 <- cor(px$expr_t2t, px$expr_patched)
lim <- range(c(px$expr_t2t, px$expr_patched))
p19a <- ggplot(px[order(g)], aes(expr_patched, expr_t2t, colour = g)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.25, colour = "grey60") +
  geom_point(size = 0.55, stroke = 0, alpha = 0.8) +
  annotate("text", x = lim[1], y = lim[2], label = sprintf("r = %.2f", r19), hjust = 0, vjust = 1, size = 6.5 / .pt,
           family = "TeX Gyre Heros") +
  scale_colour_manual(values = gcol, name = "Copies of\nthe allele") + coord_equal(xlim = lim, ylim = lim) +
  guides(colour = guide_legend(override.aes = list(size = 1.6, alpha = 1))) +
  labs(x = "Expression on patched GRCh38 (after covariates)", y = "Expression on T2T (after covariates)") + nat
txt19 <- unique(px[, .(z_indep, z_paired)])
p19b <- ggplot(px, aes(g, expr_t2t - expr_patched, colour = g)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey60") +
  geom_point(position = position_jitter(width = 0.2, height = 0, seed = 1), size = 0.45, stroke = 0, alpha = 0.7) +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = NA, linewidth = 0.3, colour = "black") +
  scale_colour_manual(values = gcol, guide = "none") +
  labs(x = "Copies of the allele", y = "Difference per person, T2T − patched") + nat
save_nat(p19a + p19b + plot_layout(widths = c(1.15, 1), guides = "collect") + plot_annotation(tag_levels = "a") &
           theme(legend.position = "right"), "fig19-paired-example", 130, 62)

# Fig 20: how many eQTLs change, the three comparisons as plain bars (replaces the stacked Fig 17 on the slide) -------
n20 <- rbindlist(lapply(names(pairs3), function(p) {
  d <- fread(file.path(GSD, paste0("changed-", p, ".tsv")))[tier %in% c("core", "main")]
  data.table(pair = p, N = nrow(d), inside = sum(d$region_class == "ffx_patched"))
}))
lab20 <- c("chm13v2.0-vs-grch38" = "T2T vs GRCh38", "grch38-ffx-vs-grch38" = "Patch alone:\npatched vs GRCh38",
           "chm13v2.0-vs-grch38-ffx" = "T2T vs patched GRCh38")
n20[, lab := factor(lab20[pair], levels = rev(lab20))]
n20[, note := fifelse(inside == 0, "none in patched regions",
              fifelse(inside == N, "all in patched regions", sprintf("%d in patched regions", inside)))]
n20[, txt := sprintf("%d   (%s)", N, note)]
n20[, col := fifelse(pair == "grch38-ffx-vs-grch38", "#E69F00", "#0072B2")]
p20 <- ggplot(n20, aes(N, lab)) + geom_col(aes(fill = col), width = 0.55) +
  lab_txt(aes(label = txt), hjust = -0.06) + scale_fill_identity() +
  scale_x_continuous(limits = c(0, max(n20$N) * 1.75), breaks = c(0, 50, 100), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = NULL) + nat +
  theme(axis.ticks = element_blank(), axis.line = element_blank(), axis.text.x = element_blank())
save_nat(p20, "fig20-changed-counts", 110, 42)

# Fig 21: the 117 changed eQTLs by segmental-duplication class (T2T SD map), with the share among all genes compared ---
gc21 <- fread(file.path(root, "results/06-compare/gene-classes.tsv"))[, .(gene_ncbi, t2t_sd_frac, t2t_sd98_frac)]
pe21 <- merge(rd("paired-effects.tsv")[, .(gene_ncbi, region_class, ch = tier %in% c("core", "main"))], gc21, by = "gene_ncbi", all.x = TRUE)
pe21[, k := fifelse(region_class == "sd", "sd98", fifelse(t2t_sd_frac >= 0.5, "sd90", "unique"))]
pe21[, k := fifelse(k == "unique", "unique", "sd")]
lv21 <- c(sd = "Segmental duplication", unique = "Unique sequence")
s21 <- pe21[, .(N = sum(ch), all = .N), by = k][, f := N / sum(N)]
s21[, lab := factor(lv21[k], levels = rev(lv21))]
s21[, txt := sprintf("%d (%.0f%%)", N, 100 * f)]
p21 <- ggplot(s21, aes(N, lab)) + geom_col(aes(fill = k), width = 0.6) +
  lab_txt(aes(label = txt), hjust = -0.1) +
  scale_fill_manual(values = c(sd = "#CC79A7", unique = "grey70"), guide = "none") +
  scale_x_continuous(limits = c(0, 95), breaks = c(0, 20, 40, 60), expand = expansion(mult = c(0, 0))) +
  labs(x = sprintf("Genes whose eQTL changes (%d in all)", s21[, sum(N)]), y = NULL) + nat +
  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank())
save_nat(p21, "fig21-sd-share", 90, 30)
fwrite(s21[, .(class = k, changed = N, all_compared = all)], file.path(dd, "changed-by-sd.tsv"), sep = "\t")

# Fig 22: reads move between copies (ASAH2 / ASAH2B), and the eQTL reverses ---------------------------------------
G22 <- "ASAH2"; P22 <- "ASAH2B"
pd22 <- dcast(rd("paralog-deltas.tsv")[pair == paste(G22, "/", P22)], sample ~ role, value.var = "delta")
r22 <- cor(pd22$gene, pd22$partner)
p22a <- ggplot(pd22, aes(gene, partner)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey80") + geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey80") +
  geom_point(size = 0.5, stroke = 0, alpha = 0.6, colour = "#CC79A7") +
  annotate("text", x = Inf, y = Inf, label = sprintf("r = %.2f", r22), hjust = 1.1, vjust = 1.5, size = 6.5 / .pt,
           family = "TeX Gyre Heros") +
  labs(x = bquote("Change in" ~ italic(.(G22)) * ", T2T − patched (CPM)"),
       y = bquote("Change in" ~ italic(.(P22)) ~ "(CPM)")) + nat
ex22 <- rd("expression-eqtl-examples.tsv")[symbol == G22]
ex22[, ref := factor(ref_lab[reference], levels = ref_lab[c("grch38-ffx", "chm13v2.0")])]
v22 <- rd("top10-changed.tsv")[symbol == G22]
ea22 <- v22$effect_allele; oa22 <- fifelse(ea22 == v22$ref_grch38, v22$alt_grch38, v22$ref_grch38)
ex22[, g := factor(round(genotype), levels = 0:2,
                   labels = c(paste0(oa22, "/", oa22), paste0(oa22, "/", ea22), paste0(ea22, "/", ea22)))]
sl22 <- unique(ex22[, .(ref, slope)])[, txt := sprintf("%s: slope %+.2f", ref, slope)]
p22b <- ggplot(ex22, aes(g, cpm, colour = ref)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.75, seed = 1), size = 0.3, stroke = 0, alpha = 0.5) +
  geom_boxplot(position = position_dodge(0.75), width = 0.6, outlier.shape = NA, fill = NA, linewidth = 0.3) +
  scale_colour_manual(values = ref_col[levels(ex22$ref)], labels = sl22[order(ref), txt], name = NULL) +
  guides(colour = guide_legend(ncol = 1)) +
  coord_cartesian(ylim = c(0, quantile(ex22$cpm, 0.995))) +
  labs(x = sprintf("Genotype at %s (effect allele %s)", v22$rsid, ea22), y = bquote(italic(.(G22)) ~ "expression (CPM)")) + nat +
  theme(legend.position = "top", legend.justification = "left")
save_nat(p22a + p22b + plot_layout(widths = c(1, 1.1)) + plot_annotation(tag_levels = "a"), "fig22-reads-move", 140, 66)

# Fig 23: top 10 changed genes as a table (R gt, rendered through grid with as_gtable, no browser needed) -------------
suppressPackageStartupMessages(library(gt))
sup <- function(x) chartr("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻", x)
pfmt <- function(l10) { e <- floor(l10); m <- 10^(l10 - e); sprintf("%.0f × 10%s", m, sup(as.character(e))) }
t10 <- rd("top10-changed.tsv")
reg23 <- c(sd = "SD", other = "Other", t2t_new_sequence = "T2T-only", ffx_patched = "FixItFelix")
comma_pos <- function(p) { x <- tstrsplit(p, ":"); sprintf("%s:%s", x[[1]], format(as.numeric(x[[2]]), big.mark = ",", trim = TRUE)) }
glob_af <- function(k, g, m) fifelse(!is.na(k), sprintf("%.2f", k),
                               fifelse(!is.na(g), sprintf("%.2f\u1d4d", g), sprintf("%.2f\u1d50", m)))
t10[, mark := paste0(fifelse(region_class == "sd", "*", ""), fifelse(sign(slope_a) != sign(slope_b), "\u2020", ""))]
tb <- t10[, .(Gene = paste0(symbol, mark), Variant = fifelse(rsid == "" | is.na(rsid), "\u2013", rsid),
              `GRCh38 position` = comma_pos(pos_grch38), `REF/ALT` = sprintf("%s/%s", ref_grch38, alt_grch38),
              Effect = effect_allele, Freq. = glob_af(eaf_1kg, eaf_gnomad, eaf_mage),
              `On T2T` = slope_a, `On patched` = slope_b, `Change (Δ)` = delta, z = z_paired,
              P = pfmt(log10p_paired),
              Replicates = fifelse(is.na(geuvadis_nominal), "n.t.", fifelse(geuvadis_nominal, "Yes", "No")))]
g23 <- gt(tb) |>
  tab_spanner("eQTL effect (SD per allele)", columns = c(`On T2T`, `On patched`, `Change (Δ)`)) |>
  tab_spanner("Paired test", columns = c(z, P)) |>
  fmt_number(columns = c(`On T2T`, `On patched`, `Change (Δ)`), decimals = 2, force_sign = TRUE) |>
  fmt_number(columns = z, decimals = 1) |>
  cols_align("right", columns = c(`GRCh38 position`, `On T2T`, `On patched`, `Change (Δ)`, z, P)) |>
  cols_align("center", columns = c(Replicates, `REF/ALT`, Effect)) |>
  cols_align("right", columns = Freq.) |>
  tab_spanner("Allele (GRCh38)", columns = c(`REF/ALT`, Effect, Freq.)) |>
  tab_style(cell_text(style = "italic"), cells_body(columns = Gene)) |>
  cols_width(Gene ~ px(62), Variant ~ px(74)) |>
  tab_source_note("* In a segmental duplication (\u2265 98% identical copies). \u2020 Reversed: the eQTL goes in opposite directions on the two references.") |>
  tab_source_note("REF/ALT as in GRCh38 and dbSNP; on T2T they are swapped at 9 of these 10 variants (CHM13 carries the other base).") |>
  tab_source_note("Effect: the allele the eQTL effects refer to. Freq.: its frequency in 1000 Genomes phase 3 (all 2,504); \u1d4d gnomAD, \u1d50 MAGE where 1000 Genomes has none.") |>
  tab_source_note("z and P: paired test (STAR). Replicates: same direction, P < 0.05 in 358 Geuvadis people not in MAGE; n.t., not testable.") |>
  tab_style(cell_text(color = "#0072B2"), cells_body(columns = `On T2T`)) |>
  tab_style(cell_text(color = "#B97A00"), cells_body(columns = `On patched`)) |>
  tab_style(cell_text(weight = "bold"), cells_body(columns = `Change (Δ)`)) |>
  tab_style(cell_text(weight = "bold"), cells_column_labels()) |>
  tab_options(table.font.names = "TeX Gyre Heros", table.font.size = px(8.5), data_row.padding = px(2.5),
              column_labels.padding = px(3), source_notes.font.size = px(7.5), table.border.top.color = "black", table.border.bottom.color = "black",
              column_labels.border.bottom.color = "black", table_body.hlines.color = "grey88",
              table.border.top.width = px(1.2), table.border.bottom.width = px(1.2), heading.border.bottom.width = px(0))
gg23 <- as_gtable(g23)
w23 <- grid::convertWidth(sum(gg23$widths), "mm", valueOnly = TRUE) + 4
h23 <- grid::convertHeight(sum(gg23$heights), "mm", valueOnly = TRUE) + 4
for (ext in c("png", "pdf")) {
  f <- file.path(fd, paste0("fig23-top10-table.", ext))
  if (ext == "png") png(f, width = w23, height = h23, units = "mm", res = 300, type = "cairo") else cairo_pdf(f, width = w23 / 25.4, height = h23 / 25.4)
  grid::grid.newpage(); grid::grid.draw(gg23); invisible(dev.off())
}
cat(sprintf("fig23 table %.0f x %.0f mm\n", w23, h23))

# Fig 24: changed genes that colocalise with a GWAS trait on one reference only (R gt) -----------------------------
cl24 <- rd("coloc-table.tsv")
cl24 <- cl24[order(on, -pmax(`PP4_chm13v2.0`, `PP4_grch38-ffx`))]
pos24 <- function(v) { x <- tstrsplit(v, ":"); sprintf("%s:%s", x[[1]], format(as.numeric(x[[2]]), big.mark = ",", trim = TRUE)) }
pf24 <- function(l10) { e <- floor(l10); m <- 10^(l10 - e); sprintf("%.0f × 10%s", m, sup(as.character(e))) }
tb24 <- cl24[, .(Gene = symbol, Trait = paste0(toupper(substr(trait, 1, 1)), substring(trait, 2)),
                 Variant = rsid, `GRCh38 position` = pos24(top), `GWAS P` = pf24(gwas_log10p),
                 `On T2T` = `PP4_chm13v2.0`, `On patched` = `PP4_grch38-ffx`,
                 `Colocalises on` = fifelse(on == "chm13v2.0", "T2T only", "Patched GRCh38 only"))]
g24 <- gt(tb24) |>
  tab_spanner("Shared variant", columns = c(Variant, `GRCh38 position`, `GWAS P`)) |>
  tab_spanner("Colocalisation (PP4)", columns = c(`On T2T`, `On patched`)) |>
  fmt_number(columns = c(`On T2T`, `On patched`), decimals = 2) |>
  cols_align("right", columns = c(`GRCh38 position`, `GWAS P`, `On T2T`, `On patched`)) |>
  tab_style(cell_text(style = "italic"), cells_body(columns = Gene)) |>
  tab_style(cell_text(weight = "bold", color = "#0072B2"), cells_body(columns = `On T2T`, rows = `On T2T` >= 0.8)) |>
  tab_style(cell_text(weight = "bold", color = "#B97A00"), cells_body(columns = `On patched`, rows = `On patched` >= 0.8)) |>
  tab_style(cell_text(weight = "bold"), cells_column_labels()) |>
  cols_width(Gene ~ px(92)) |>
  tab_source_note("Shared variant: the variant with the most combined eQTL and GWAS evidence, on the reference where the gene colocalises.") |>
  tab_source_note("PP4: probability that the eQTL and the GWAS signal share one causal variant (coloc); ≥ 0.8 in bold.") |>
  tab_options(table.font.names = "TeX Gyre Heros", table.font.size = px(9), data_row.padding = px(3),
              column_labels.padding = px(3), source_notes.font.size = px(7.5), table.border.top.color = "black",
              table.border.bottom.color = "black", column_labels.border.bottom.color = "black",
              table_body.hlines.color = "grey88", table.border.top.width = px(1.2), table.border.bottom.width = px(1.2))
gg24 <- as_gtable(g24)
w24 <- grid::convertWidth(sum(gg24$widths), "mm", valueOnly = TRUE) + 4
h24 <- grid::convertHeight(sum(gg24$heights), "mm", valueOnly = TRUE) + 4
for (ext in c("png", "pdf")) {
  f <- file.path(fd, paste0("fig24-gwas-table.", ext))
  if (ext == "png") png(f, width = w24, height = h24, units = "mm", res = 300, type = "cairo") else cairo_pdf(f, width = w24 / 25.4, height = h24 / 25.4)
  grid::grid.newpage(); grid::grid.draw(gg24); invisible(dev.off())
}
cat(sprintf("fig24 table %.0f x %.0f mm\n", w24, h24))

# Fig 25: RNA coverage along a paralog pair on both references, 0 vs 2 copies of the eQTL allele (reads moving) -------
cv <- rd("coverage-example.tsv")[genotype %in% c(0, 2)]
cv[, ref := factor(ref_lab[ref], levels = ref_lab[c("grch38-ffx", "chm13v2.0")])]
cv[, g := factor(genotype, levels = c(0, 2), labels = c("0 copies", "2 copies"))]
genes25 <- unique(cv$gene)
v25 <- rd("top10-changed.tsv")[symbol == genes25[1]]
ea25 <- v25$effect_allele; oa25 <- fifelse(v25$effect_allele == v25$ref_grch38, v25$alt_grch38, v25$ref_grch38)
rs25 <- v25$rsid; gt0 <- paste0(oa25, "/", oa25); gt2 <- paste0(ea25, "/", ea25)
n25 <- unique(rd("coverage-example.tsv")[, .(genotype, n)])
rat25 <- cv[, .(tot = sum(cpm)), by = .(ref, gene, genotype)][, .(r = tot[genotype == 2] / tot[genotype == 0]), by = .(ref, gene)]
cov_col <- function(gn, ylab = NULL) {
  d <- cv[gene == gn]
  rr <- rat25[gene == gn]
  rr[, txt := sprintf("%s vs %s: %.2f\u00d7", gt2, gt0, r)]
  ggplot(d, aes(pos_kb, cpm)) +
    geom_area(data = d[g == "0 copies"], fill = "grey75", alpha = 0.8) +
    geom_area(data = d[g == "2 copies"], fill = "#0072B2", alpha = 0.45) +
    geom_line(aes(colour = g), linewidth = 0.3) +
    geom_text(data = rr, aes(x = Inf, y = Inf, label = txt), hjust = 1.05, vjust = 1.4, size = 6.5 / .pt,
              family = "TeX Gyre Heros", fontface = fifelse(abs(log(rr$r)) > log(1.3), "bold", "plain"), inherit.aes = FALSE) +
    facet_grid(ref ~ ., switch = "y") +
    scale_colour_manual(values = c("0 copies" = "grey45", "2 copies" = "#0072B2"),
                        labels = sprintf("%s (%d people)", c(gt0, gt2), n25[match(c(0, 2), genotype), n]),
                        name = sprintf("Genotype at %s (effect allele %s):", rs25, ea25)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.15))) +
    labs(title = gn, x = "Exons joined, 5′ to 3′ (kb)", y = ylab) + nat +
    theme(plot.title = element_text(size = 7, face = "bold.italic", hjust = 0.5), strip.placement = "outside",
          strip.text.y.left = element_text(angle = 90, face = "bold"), panel.spacing = unit(5, "pt"))
}
p25 <- (cov_col(genes25[1], "RNA coverage (reads per million mapped)") + cov_col(genes25[2]) +
  plot_layout(guides = "collect")) & theme(legend.position = "bottom")
save_nat(p25, "fig25-coverage", 150, 72)

# Fig 25s (supplementary): the same coverage with all three genotypes, heterozygotes included -------------------------
cs <- rd("coverage-example.tsv")
cs[, ref := factor(ref_lab[ref], levels = ref_lab[c("grch38-ffx", "chm13v2.0")])]
vs <- rd("top10-changed.tsv")[symbol == unique(cs$gene)[1]]
eas <- vs$effect_allele; oas <- fifelse(eas == vs$ref_grch38, vs$alt_grch38, vs$ref_grch38)
gl <- c(paste0(oas, "/", oas), paste0(oas, "/", eas), paste0(eas, "/", eas))
ns <- unique(cs[, .(genotype, n)])[order(genotype)]
cs[, g := factor(genotype, levels = 0:2, labels = sprintf("%s (%d people)", gl, ns$n))]
rs <- cs[, .(tot = sum(cpm)), by = .(ref, gene, genotype)]
rs <- rs[, .(txt = sprintf("%s : %s : %s = 1 : %.2f : %.2f", gl[1], gl[2], gl[3], tot[genotype == 1] / tot[genotype == 0],
                         tot[genotype == 2] / tot[genotype == 0])), by = .(ref, gene)]
p25s <- ggplot(cs, aes(pos_kb, cpm, colour = g)) +
  geom_line(linewidth = 0.35) +
  geom_text(data = rs, aes(x = Inf, y = Inf, label = txt), hjust = 1.03, vjust = 1.4, size = 6 / .pt,
            family = "TeX Gyre Heros", inherit.aes = FALSE) +
  facet_grid(ref ~ gene, scales = "free", switch = "y") +
  scale_colour_manual(values = setNames(c("grey55", "#56B4E9", "#0072B2"), levels(cs$g)),
                      name = sprintf("Genotype at %s (effect allele %s):", vs$rsid, eas)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Exons joined, 5′ to 3′ (kb)", y = "RNA coverage (reads per million mapped)") + nat +
  theme(strip.text.x = element_text(face = "bold.italic"), strip.placement = "outside",
        panel.spacing = unit(6, "pt"), legend.position = "bottom")
save_nat(p25s, "fig25s-coverage-all-genotypes", 150, 76)
