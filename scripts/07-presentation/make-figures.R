# Presentation figures 2-10 (ggplot2) from results/07-presentation/data/ (prepare-figure-data.py).
#
#   pixi run Rscript scripts/07-presentation/make-figures.R
#
# Writes results/07-presentation/figures/figNN-*.png (300 dpi). Figure 1 (pipeline) is Graphviz: pipeline.dot.
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

# Fig 2: cohorts ---------------------------------------------------------------------------------------------
co <- unique(rd("cohorts.tsv")[!is.na(population) & population != ""], by = c("cohort", "sample"))
co[, cohort := factor(ifelse(cohort == "mage", "MAGE (discovery), n = 731", "Geuvadis,\nn = 451"),
                      levels = c("MAGE (discovery), n = 731", "Geuvadis,\nn = 451"))]
pc <- co[, .N, by = .(cohort, population, superpop)]
pc[, population := reorder(population, as.integer(factor(superpop)) * 1000 + N)]
p2a <- ggplot(pc, aes(population, N, fill = superpop)) + geom_col(width = 0.75) +
  facet_grid(. ~ cohort, scales = "free_x", space = "free_x") + scale_fill_manual(values = sp_col, name = "Super-population") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Individuals", title = "Two lymphoblastoid cohorts with 1000 Genomes genotypes") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 9))
rp <- rd("cohorts.tsv")
rp[, cohort := ifelse(cohort == "mage", "MAGE: 2x150 bp", "Geuvadis: 2x76 bp")]
p2b <- ggplot(rp, aes(read_pairs / 1e6, fill = cohort)) + geom_histogram(bins = 50, alpha = 0.8, position = "identity") +
  scale_fill_manual(values = c("MAGE: 2x150 bp" = "#0072B2", "Geuvadis: 2x76 bp" = "#56B4E9"), name = NULL) +
  labs(x = "Read pairs per run (millions)", y = "Runs", title = "Sequencing depth per run") + coord_cartesian(xlim = c(0, 100))
save(p2a / p2b + plot_layout(heights = c(1.4, 1)), "fig02-cohorts.png", 11, 8.5)

# Fig 3: mapping QC ------------------------------------------------------------------------------------------
qc <- rd("mapping-qc.tsv")
qm <- melt(qc, id.vars = c("run", "sample"), measure.vars = patterns(unique = "^unique_pct_", multi = "^multi_pct_"),
           variable.name = "ref")
qm[, ref := factor(ref_lab[c("chm13v2.0", "grch38", "grch38-ffx")[as.integer(ref)]], levels = ref_lab)]
qm2 <- melt(qm, id.vars = c("run", "ref"), measure.vars = c("unique", "multi"))
qm2[, variable := factor(ifelse(variable == "unique", "Uniquely mapped reads (%)", "Multi-mapped reads (%)"),
                         levels = c("Uniquely mapped reads (%)", "Multi-mapped reads (%)"))]
med <- qm2[, .(m = median(value)), by = .(ref, variable)]
p3 <- ggplot(qm2, aes(ref, value, colour = ref)) +
  geom_line(aes(group = run), colour = "grey85", linewidth = 0.2) + geom_boxplot(outlier.size = 0.4, width = 0.5, fill = NA) +
  geom_text(data = med, aes(y = m, label = sprintf("%.2f", m)), nudge_x = 0.38, size = 3.6, colour = "black") +
  facet_wrap(~variable, scales = "free_y") + scale_colour_manual(values = ref_col, guide = "none") +
  labs(x = NULL, y = NULL, title = "Same reads, three references: T2T places more reads uniquely",
       subtitle = "MAGE, 779 runs; grey lines connect the same run; numbers are medians")
save(p3, "fig03-mapping-qc.png", 10, 5)

# Fig 4: variant matching ------------------------------------------------------------------------------------
vm <- rd("variant-matching.tsv")[, .(n = sum(n)), by = .(side, category, type)]
vm[, label := fcase(category == "matched", "Matched (same variant on both)",
                    category == "matched_maf_discordant", "Paired but MAF differs (paralogous)",
                    category == "lifted_not_called_grch38", "T2T variant, not called on GRCh38",
                    category == "no_liftover", "In T2T-only sequence",
                    category == "not_called_t2t", "GRCh38 variant, not called on T2T")]
vm[, label := factor(label, levels = rev(c("Matched (same variant on both)", "Paired but MAF differs (paralogous)",
                                           "T2T variant, not called on GRCh38", "In T2T-only sequence",
                                           "GRCh38 variant, not called on T2T")))]
p4 <- ggplot(vm, aes(n / 1e6, label, fill = type)) + geom_col(width = 0.7) +
  geom_text(data = vm[, .(n = sum(n)), by = label], aes(n / 1e6, label, label = sprintf("%.2f M", n / 1e6)),
            inherit.aes = FALSE, hjust = -0.1, size = 3.8) +
  scale_fill_manual(values = c(snp = "#0072B2", indel = "#56B4E9"), name = NULL, labels = c(indel = "Indel", snp = "SNP")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Variants with MAF >= 1% in the 1,089 expression individuals (millions)", y = NULL,
       title = "Pairing variants between the T2T and GRCh38 genotype panels",
       subtitle = "The main comparison uses only matched variants,\nso differences come from expression, not genotypes")
save(p4, "fig04-variant-matching.png", 11, 4.8)

# Fig 5: eGenes ----------------------------------------------------------------------------------------------
eg <- rd("egenes.tsv")
eg[, ref := factor(ref_lab[reference], levels = ref_lab)]
eg[, method := factor(ifelse(method == "star", "STAR counts", "salmon counts"), levels = c("STAR counts", "salmon counts"))]
p5 <- ggplot(eg, aes(ref, egenes, fill = ref)) + geom_col(width = 0.65) +
  geom_text(aes(label = comma(egenes)), vjust = -0.4, size = 4) +
  geom_text(aes(y = 600, label = paste0("of ", comma(genes_tested), "\ntested")), colour = "white", size = 3.2, lineheight = 0.9) +
  facet_wrap(~method) + scale_fill_manual(values = ref_col, guide = "none") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "eGenes (FDR 5%)", title = "About 15,000 eGenes on every reference",
       subtitle = "MAGE, 731 individuals, chr1-22, matched variants\n(the MAGE paper reports ~15,000 on GRCh38)")
save(p5, "fig05-egenes.png", 9, 5)

# Fig 6: positive control ------------------------------------------------------------------------------------
pcq <- rd("positive-control.tsv")
pcq[, ref := factor(ref_lab[reference], levels = ref_lab)]
pcq[, symbol := factor(symbol, levels = pcq[reference == "grch38"][order(qval)]$symbol)]
p6 <- ggplot(pcq, aes(ref, -log10(qval), group = symbol, colour = symbol)) + geom_line(linewidth = 1) + geom_point(size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey50") +
  geom_text_repel(data = pcq[reference == "grch38"], aes(label = symbol), nudge_x = -0.3, direction = "y", size = 4, show.legend = FALSE) +
  scale_y_log10() + scale_colour_brewer(palette = "Dark2", guide = "none") +
  labs(x = NULL, y = "-log10 q (log scale)", title = "Positive control: FixItFelix's chr21 artefacts reappear",
       subtitle = "Inflated on unpatched GRCh38; T2T agrees with the patched reference.\nDashed line: q = 0.05")
save(p6, "fig06-positive-control.png", 8, 5.5)

# Fig 7: genome-wide same-variant comparison -----------------------------------------------------------------
ef <- rd("effects.tsv")
ef[, cls := factor(cls_lab[region_class], levels = cls_lab)]
ef[, grp := fifelse(robust, "Robust change (STAR and salmon)", fifelse(changed, "Changed with STAR only", "No change"))]
lab <- ef[robust == TRUE][order(p_diff)][1:18]
p7 <- ggplot(ef[order(grp == "No change", decreasing = TRUE)], aes(slope_b, slope_a)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60") + geom_hline(yintercept = 0, colour = "grey85") +
  geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(data = ef[grp == "No change"], colour = "grey75", size = 0.5, alpha = 0.5) +
  geom_point(data = ef[grp == "Changed with STAR only"], colour = "grey40", size = 1.2, shape = 1) +
  geom_point(data = ef[robust == TRUE], aes(colour = cls), size = 2.4) +
  geom_text_repel(data = lab, aes(label = symbol, colour = cls), size = 3.3, max.overlaps = 30, show.legend = FALSE) +
  scale_colour_manual(values = cls_col, name = NULL, drop = TRUE) +
  labs(x = "Effect size on patched GRCh38 (slope, same variant)", y = "Effect size on T2T-CHM13 (slope)",
       title = "Same gene, same variant: where do the two references disagree?",
       subtitle = sprintf("%s genes. Grey dots: no change; open circles: changed with STAR only;\ncoloured: %d robust changes (STAR and salmon), none in FixItFelix regions",
                          comma(nrow(ef)), sum(ef$robust)))
save(p7, "fig07-effect-comparison.png", 9.5, 8.5)

# Fig 8: the robust genes ------------------------------------------------------------------------------------
rg <- rd("robust-genes.tsv")
rg[, cls := factor(cls_lab[region_class], levels = cls_lab)]
rg[, symbol := reorder(symbol, z_diff)]
rg[, side := ifelse(stronger_on == "chm13v2.0", "Stronger on T2T", "Stronger on patched GRCh38")]
p8a <- ggplot(rg, aes(z_diff, symbol, colour = cls)) + geom_vline(xintercept = 0, colour = "grey70") +
  geom_segment(aes(x = 0, xend = z_diff, yend = symbol), linewidth = 0.6) + geom_point(size = 2.2) +
  scale_colour_manual(values = cls_col, name = NULL, drop = TRUE) +
  labs(x = "Effect difference z (T2T - patched GRCh38), STAR", y = NULL, title = "57 genes changed beyond the patch") +
  theme(axis.text.y = element_text(size = 7.5))
cnt <- rg[, .N, by = .(cls, side)]
p8b <- ggplot(cnt, aes(side, N, fill = cls)) + geom_col(width = 0.6) + geom_text(aes(label = N), position = position_stack(vjust = 0.5), colour = "white") +
  scale_fill_manual(values = cls_col, guide = "none") + labs(x = NULL, y = "Genes", title = "By region class") +
  theme(axis.text.x = element_text(size = 10))
rg[, partner_ok := !is.na(partner_r) & partner_r < -0.5]
p8c <- ggplot(rg, aes(partner_r, fill = partner_ok)) + geom_histogram(binwidth = 0.1, boundary = 0) +
  scale_fill_manual(values = c(`TRUE` = "#CC79A7", `FALSE` = "grey70"), labels = c(`TRUE` = "clear partner (r < -0.5)", `FALSE` = "no single partner"), name = NULL) +
  labs(x = "Correlation with the gene that takes its reads", y = "Genes", title = "Read destination found for 33/57") +
  theme(legend.position = "bottom")
save(p8a | (p8b / p8c), "fig08-robust-genes.png", 12, 10)

# Fig 9: paralog swaps ---------------------------------------------------------------------------------------
pd <- rd("paralog-deltas.tsv")
pw <- dcast(pd, pair + sample ~ role, value.var = "delta")
r9 <- pw[, .(r = cor(gene, partner)), by = pair]
pw <- merge(pw, r9, by = "pair")
pw[, pair := factor(paste0(pair, "   (r = ", sprintf("%.2f", r), ")"))]
p9 <- ggplot(pw, aes(gene, partner)) + geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_abline(slope = -1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(size = 0.8, alpha = 0.5, colour = "#CC79A7") + geom_smooth(method = "lm", formula = y ~ x, colour = "black", linewidth = 0.6, se = FALSE) +
  facet_wrap(~pair, scales = "free") +
  labs(x = "Change in gene (T2T - patched GRCh38, CPM)", y = "Change in partner (CPM)",
       title = "Where the reads go: one gene loses what its partner gains",
       subtitle = "Each dot is one MAGE individual; dashed line: reads exactly swapped.\nLOC128966601 lies in T2T-only sequence; RPSA2 is a processed pseudogene of RPSA")
save(p9, "fig09-paralog-swaps.png", 10, 8)

# Fig 10: eQTL examples --------------------------------------------------------------------------------------
ex <- rd("eqtl-examples.tsv")
ex[, ref := factor(ref_lab[reference], levels = ref_lab)]
ex[, g := factor(round(genotype), levels = 0:2, labels = c("0", "1", "2"))]
ex[, panel := factor(fcase(symbol == "PWP2", "PWP2: FixItFelix artefact (GRCh38 vs patched)",
                           symbol == "CCZ1", "CCZ1: effect reverses (paralog CCZ1B)",
                           symbol == "SMN2", "SMN2: stronger on T2T (SMN1/SMN2)"),
                     levels = c("PWP2: FixItFelix artefact (GRCh38 vs patched)", "CCZ1: effect reverses (paralog CCZ1B)",
                                "SMN2: stronger on T2T (SMN1/SMN2)"))]
p10 <- ggplot(ex[!is.na(g)], aes(g, expression, fill = ref)) +
  geom_boxplot(outlier.size = 0.3, width = 0.7, position = position_dodge(0.8)) +
  facet_wrap(~panel, scales = "free_y") + scale_fill_manual(values = ref_col, name = NULL) +
  labs(x = "Genotype at the lead variant (copies of the alternative allele)", y = "Normalised expression",
       title = "The same people, the same variant, a different answer",
       subtitle = "MAGE, 731 individuals; expression from each reference's alignments")
save(p10, "fig10-eqtl-examples.png", 13, 5.2)
cat("figures written to", fd, "\n")
