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
            color_value = pmin(-log10(FDR), 6))

make_plot <- function(ranks, full = FALSE) {
 z <- long |> filter(display_rank %in% ranks)
 names_in_order <- d$display_label[d$display_rank %in% ranks]
 z$display_label <- factor(z$display_label, levels = rev(names_in_order))
 # Plotmath bold labels emphasize the two discussed pathways without changing their rank.
 axis_labels <- lapply(rev(names_in_order), function(s) {
   if (s %in% c('Proteasome', 'Mitophagy')) bquote(bold(.(s))) else s
 })
 names(axis_labels) <- rev(names_in_order)
 axis_labels <- as.expression(axis_labels)
 ggplot(z, aes(NES, display_label)) +
   geom_point(aes(color = color_value), size = if (full) 2.7 else 3.15) +
   facet_grid(. ~ Dataset) +
   scale_y_discrete(labels = axis_labels, expand = expansion(add = .65)) +
   scale_x_continuous(limits = c(1.3, 2.9), breaks = c(1.5, 2, 2.5), expand = expansion(mult = 0)) +
   scale_color_gradientn(colors = c('#91A8CB', '#788BB7', '#5E427D', '#AE355E', '#C84A35'),
     limits = c(-log10(.05), 6), breaks = c(-log10(.05), 2, 3, 6),
     labels = c('0.05', '0.01', '0.001', '\u226410\u207b\u2076'), name = 'FDR',
     guide = guide_colorbar(direction = 'horizontal', title.position = 'left',
       barwidth = grid::unit(65, 'mm'), barheight = grid::unit(2.5, 'mm'), ticks = TRUE)) +
   labs(x = 'Normalized enrichment score (NES)', y = NULL) +
   theme_classic(base_family = 'Arial', base_size = 10) +
   theme(axis.text.y = element_text(size = if (full) 9.4 else 11, color = '#202020', margin = margin(r = 8)),
         axis.text.x = element_text(size = 9, color = '#303030'),
         axis.title.x = element_text(size = 10, margin = margin(t = 7)),
         axis.line.y = element_blank(), axis.ticks.y = element_blank(),
         axis.line.x = element_line(linewidth = .35), axis.ticks.x = element_line(linewidth = .3),
         panel.grid.major.y = element_line(color = '#EEEEEE', linewidth = .28),
         panel.grid.major.x = element_blank(), strip.background = element_blank(),
         strip.text = element_text(size = 11, face = 'bold', margin = margin(b = 9)),
         panel.spacing.x = grid::unit(7, 'mm'), legend.position = 'bottom',
         legend.margin = margin(t = 3), legend.text = element_text(size = 8.5),
         legend.title = element_text(size = 9), plot.margin = margin(8, 8, 5, 5))
}
save_plot <- function(p, name, w, h) {
 ggsave(file.path(out, paste0(name, '.png')), p, width = w, height = h, dpi = 320, bg = 'white', device = ragg::agg_png)
 ggsave(file.path(out, paste0(name, '.pdf')), p, width = w, height = h, bg = 'white', device = cairo_pdf)
}
save_plot(make_plot(1:40, TRUE), 'Fig06_shared_pathways_dotplot', 8.1, 10.4)
save_plot(make_plot(1:20), 'Fig06B_shared_pathways_ranks01_20', 10.9, 5.8)
save_plot(make_plot(21:40), 'Fig06C_shared_pathways_ranks21_40', 10.9, 5.8)
message('Verified 40 pathways / 80 unchanged NES-FDR pairs. Proteasome rank 14; Mitophagy rank 28.')
