options(stringsAsFactors = FALSE)
.libPaths(unique(c(Sys.getenv('GSE_REPRO_R_LIB', 'C:/Users/User/Documents/R/win-library/4.4'), .libPaths())))
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
args <- commandArgs(TRUE)
if (length(args) != 2L) stop('Usage: Rscript 26_shared_pathway_dotplot.R Table07.csv figure_directory')
d <- read.csv(args[1], check.names = FALSE)
out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
stopifnot(nrow(d) == 40L, !anyDuplicated(d$key), all(d$analysis == 'within15'),
          all(d$sc_FDR < .05), all(d$bulk_FDR < .05), all(d$sc_NES > 0), all(d$bulk_NES > 0))
# Display order only. The larger of the two FDRs is NOT a combined P value.
d <- d |> arrange(pmax(sc_FDR, bulk_FDR), Pathway) |> mutate(display_rank = row_number())
stopifnot(d$display_rank[d$key == 'proteasome'] == 14L,
          d$display_rank[d$key == 'mitophagy'] == 28L)
label_map <- c(
 'Human immunodeficiency virus 1 infection' = 'HIV-1 infection',
 'Human T-cell leukemia virus 1 infection' = 'HTLV-1 infection',
 'Kaposi sarcoma-associated herpesvirus infection' = 'KSHV infection',
 'Human cytomegalovirus infection' = 'Cytomegalovirus infection'
)
d$display_label <- d$Pathway
ii <- match(d$Pathway, names(label_map))
d$display_label[!is.na(ii)] <- unname(label_map[ii[!is.na(ii)]])
d$display_label <- sub(' signaling pathway$', ' signaling', d$display_label)
long <- bind_rows(
 d |> transmute(key, Pathway, display_label, display_rank, Dataset = 'Relapse scRNA-seq', NES = sc_NES, FDR = sc_FDR),
 d |> transmute(key, Pathway, display_label, display_rank, Dataset = 'STAT5B N642H bulk', NES = bulk_NES, FDR = bulk_FDR)
) |> mutate(Dataset = factor(Dataset, levels = c('Relapse scRNA-seq', 'STAT5B N642H bulk')),
            FDR_label = ifelse(FDR < .001, formatC(FDR, format = 'e', digits = 2),
                               formatC(FDR, format = 'f', digits = 4)),
            stars = ifelse(FDR < .001, '***', ifelse(FDR < .01, '**', '*')))

make_plot <- function(ranks, full = FALSE) {
 z <- long |> filter(display_rank %in% ranks)
 names_in_order <- d$display_label[d$display_rank %in% ranks]
 z$row <- length(ranks) + 1 - match(z$display_rank, ranks)
 # Plotmath bold labels emphasize the two discussed pathways without changing their rank.
 axis_labels <- lapply(rev(names_in_order), function(s) {
   if (s %in% c('Proteasome', 'Mitophagy')) bquote(bold(.(s))) else s
 })
 names(axis_labels) <- rev(names_in_order)
 axis_labels <- as.expression(axis_labels)
 n <- length(ranks)
 text_size <- if (full) 3.0 else 3.25
 ggplot(z, aes(NES, row)) +
   geom_segment(aes(x = 1.3, xend = 3.18, yend = row), color = '#EAEAEA', linewidth = .28) +
   geom_point(aes(color = Dataset), size = 2.65) +
   geom_text(aes(label = sprintf('%.2f', NES)), nudge_x = .085, hjust = 0,
             size = text_size, color = '#202020') +
   geom_text(aes(x = 4.12, label = FDR_label), hjust = 1,
             size = text_size, color = '#202020') +
   geom_text(aes(x = 4.21, label = stars), hjust = 0, size = text_size,
             fontface = 'bold', color = '#202020') +
   annotate('segment', x = 3.39, xend = 3.39, y = .45, yend = n + .48,
            color = '#C9CDD2', linewidth = .35) +
   annotate('text', x = 4.0, y = n + 1.03, label = 'FDR', fontface = 'bold', size = 3.55) +
   annotate('segment', x = 1.3, xend = 3.18, y = .3, yend = .3, linewidth = .35) +
   geom_text(data = data.frame(x = c(1.5, 2, 2.5, 3), y = -.25, label = c('1.5', '2.0', '2.5', '3.0')),
             aes(x, y, label = label), inherit.aes = FALSE, size = 3) +
   annotate('text', x = 2.25, y = -1.15, label = 'NES', size = 3.6) +
   facet_grid(. ~ Dataset) +
   scale_y_continuous(breaks = seq_len(n), labels = axis_labels, limits = c(-1.6, n + 1.6), expand = expansion(mult = 0)) +
   scale_x_continuous(limits = c(1.27, 4.62), breaks = NULL, expand = expansion(mult = 0)) +
   scale_color_manual(values = c('Relapse scRNA-seq' = '#0072B2', 'STAT5B N642H bulk' = '#D55E00'), guide = 'none') +
   labs(x = NULL, y = NULL,
        caption = 'NES, normalized enrichment score.  * FDR < 0.05; ** FDR < 0.01; *** FDR < 0.001') +
   theme_classic(base_family = 'Arial', base_size = 10) +
   theme(axis.text.y = element_text(size = if (full) 9.3 else 10.7, color = '#202020', margin = margin(r = 8)),
         axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
         panel.grid = element_blank(), strip.background = element_blank(),
         strip.text = element_text(size = 11, face = 'bold', margin = margin(b = 9)),
         panel.spacing.x = grid::unit(6, 'mm'),
         plot.caption = element_text(size = 8.5, hjust = 0, margin = margin(t = 5)),
         plot.margin = margin(8, 8, 5, 5))
}
save_plot <- function(p, name, w, h) {
 built <- ggplot_build(p)
 stopifnot(nrow(built$data[[2]]) == nrow(p$data),
           isTRUE(all.equal(built$data[[2]]$x, p$data$NES, tolerance = 0)),
           identical(built$data[[3]]$label, sprintf('%.2f', p$data$NES)),
           identical(built$data[[4]]$label, p$data$FDR_label),
           identical(built$data[[5]]$label, p$data$stars))
 ggsave(file.path(out, paste0(name, '.png')), p, width = w, height = h, dpi = 320, bg = 'white', device = ragg::agg_png)
 ggsave(file.path(out, paste0(name, '.pdf')), p, width = w, height = h, bg = 'white', device = cairo_pdf)
}
save_plot(make_plot(1:40, TRUE), 'Fig06_shared_pathways_dotplot', 10.8, 10.4)
save_plot(make_plot(1:20), 'Fig06B_shared_pathways_ranks01_20', 11.5, 5.8)
save_plot(make_plot(21:40), 'Fig06C_shared_pathways_ranks21_40', 11.5, 5.8)
message('Verified 40 pathways / 80 unchanged NES-FDR pairs. Proteasome rank 14; Mitophagy rank 28.')
