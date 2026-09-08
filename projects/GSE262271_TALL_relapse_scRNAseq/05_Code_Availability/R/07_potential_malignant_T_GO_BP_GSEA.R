options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(args)) args[[1]] else "outputs/GSE262271_TALL_relapse_scRNAseq_final"
output_root <- normalizePath(output_root, winslash = "/", mustWork = TRUE)

user_lib <- Sys.getenv("R_LIBS_USER", unset = "C:/Users/User/Documents/R/win-library/4.4")
if (dir.exists(user_lib)) .libPaths(unique(c(normalizePath(user_lib, winslash = "/"), .libPaths())))
required <- c("dplyr", "ggplot2", "fgsea", "AnnotationDbi", "org.Hs.eg.db", "GO.db")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "), call. = FALSE)

csv_dir <- file.path(output_root, "csv", "potential_malignant_T_cells")
png_dir <- file.path(output_root, "png", "potential_malignant_T_cells")
deg_file <- file.path(csv_dir, "pseudobulk_DEG_relapse_vs_diagnosis_all_genes.csv")
deg <- utils::read.csv(deg_file, check.names = FALSE)
deg <- deg[is.finite(deg$rank_stat) & !is.na(deg$gene) & nzchar(deg$gene), , drop = FALSE]
deg <- deg[order(-abs(deg$rank_stat)), , drop = FALSE]
deg <- deg[!duplicated(deg$gene), , drop = FALSE]
rank_stats <- deg$rank_stat
names(rank_stats) <- deg$gene
rank_stats <- sort(rank_stats, decreasing = TRUE)

valid_symbols <- intersect(names(rank_stats), AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "SYMBOL"))
go_map <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = valid_symbols,
  columns = c("GOALL", "ONTOLOGYALL"),
  keytype = "SYMBOL"
))
go_map <- go_map |>
  dplyr::filter(!is.na(GOALL), ONTOLOGYALL == "BP") |>
  dplyr::distinct(SYMBOL, GOALL)
go_terms <- suppressMessages(AnnotationDbi::select(
  GO.db::GO.db,
  keys = unique(go_map$GOALL),
  columns = c("TERM", "ONTOLOGY"),
  keytype = "GOID"
))
go_terms <- go_terms |>
  dplyr::filter(ONTOLOGY == "BP") |>
  dplyr::distinct(GOID, TERM)
go_joined <- go_map |>
  dplyr::inner_join(go_terms, by = c("GOALL" = "GOID"))

go_sets_df <- go_joined |>
  dplyr::group_by(GOALL, TERM) |>
  dplyr::summarise(genes = list(sort(unique(SYMBOL))), .groups = "drop") |>
  dplyr::mutate(
    pathway = paste0(GOALL, " ", TERM),
    ranked_genes_n = lengths(genes),
    tested = ranked_genes_n >= 10 & ranked_genes_n <= 500
  )
tested_df <- go_sets_df |>
  dplyr::filter(tested)
pathways <- stats::setNames(tested_df$genes, tested_df$pathway)

set.seed(20260716)
fg <- fgsea::fgseaMultilevel(
  pathways = pathways,
  stats = rank_stats,
  minSize = 10,
  maxSize = 500,
  eps = 0
)
fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
fg <- as.data.frame(fg)

audit <- go_sets_df |>
  dplyr::transmute(
    gene_set_source = "GO_BP",
    gene_set_id = GOALL,
    pathway_name = TERM,
    pathway,
    ranked_genes_n,
    tested,
    untested_reason = dplyr::case_when(
      tested ~ "",
      ranked_genes_n < 10 ~ "fewer_than_10_ranked_genes",
      ranked_genes_n > 500 ~ "more_than_500_ranked_genes",
      TRUE ~ "not_tested"
    ),
    linked_genes = vapply(genes, paste, collapse = ";", FUN.VALUE = character(1))
  )

result <- audit |>
  dplyr::left_join(fg, by = "pathway") |>
  dplyr::mutate(
    analysis_population = "potential malignant T cells (CNV-high top 25%)",
    comparison = "relapse_vs_diagnosis_paired_pseudobulk_N3",
    direction = dplyr::case_when(
      !tested ~ "untested",
      NES > 0 ~ "relapse_up",
      NES < 0 ~ "relapse_down",
      TRUE ~ "flat"
    ),
    significant_FDR_0.05 = tested & !is.na(padj) & padj < 0.05
  ) |>
  dplyr::select(
    analysis_population, comparison, gene_set_source, gene_set_id, pathway_name,
    pathway, ranked_genes_n, tested, untested_reason, pval, padj, log2err, ES,
    NES, size, direction, significant_FDR_0.05, leadingEdge, linked_genes
  ) |>
  dplyr::arrange(dplyr::desc(tested), padj, dplyr::desc(abs(NES)))

utils::write.csv(
  result,
  file.path(csv_dir, "GO_BP_GSEA_all_terms_relapse_vs_diagnosis.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

sig <- result |>
  dplyr::filter(significant_FDR_0.05) |>
  dplyr::arrange(padj, dplyr::desc(abs(NES)))
utils::write.csv(
  sig,
  file.path(csv_dir, "GO_BP_GSEA_significant_terms_FDR0.05.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

summary <- data.frame(
  analysis_population = "potential malignant T cells (CNV-high top 25%)",
  GO_BP_terms_total = nrow(result),
  GO_BP_terms_tested = sum(result$tested),
  significant_FDR_0.05_total = nrow(sig),
  significant_relapse_up = sum(sig$direction == "relapse_up"),
  significant_relapse_down = sum(sig$direction == "relapse_down"),
  stringsAsFactors = FALSE
)
utils::write.csv(
  summary,
  file.path(csv_dir, "GO_BP_GSEA_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

top_each_direction <- function(x, n = 20) {
  if (!nrow(x)) return(x)
  head(x[order(x$padj, -abs(x$NES)), , drop = FALSE], n)
}
plot_df <- rbind(
  top_each_direction(sig[sig$direction == "relapse_up", , drop = FALSE], 20),
  top_each_direction(sig[sig$direction == "relapse_down", , drop = FALSE], 20)
)
if (!nrow(plot_df)) stop("No significant GO BP terms to plot.", call. = FALSE)
plot_df$minus_log10_FDR <- -log10(plot_df$padj)
plot_df$display_label <- paste0(
  plot_df$pathway_name,
  ifelse(plot_df$direction == "relapse_up", " [Up]", " [Down]")
)
plot_df$display_label <- factor(
  plot_df$display_label,
  levels = plot_df$display_label[order(plot_df$NES, decreasing = FALSE)]
)

p <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(x = NES, y = display_label, size = size, color = minus_log10_FDR)
) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.35, color = "grey55") +
  ggplot2::geom_point(alpha = 0.9) +
  ggplot2::scale_color_gradient(low = "#4C78A8", high = "#B2433F") +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.10, 0.10))) +
  ggplot2::labs(
    title = "Significant GO Biological Process GSEA in potential malignant T cells",
    subtitle = "Relapse vs diagnosis; paired pseudobulk (N = 3); top 20 per direction by FDR",
    x = "Normalized enrichment score (NES)",
    y = NULL,
    size = "Gene-set size",
    color = expression(-log[10]("FDR"))
  ) +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(
    axis.text = ggplot2::element_text(color = "black"),
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "right",
    panel.grid.major.x = ggplot2::element_line(color = "grey90", linewidth = 0.3)
  )
ggplot2::ggsave(
  file.path(png_dir, "GO_BP_GSEA_significant_terms_dotplot.png"),
  p,
  width = 10.4,
  height = max(7.0, 0.31 * nrow(plot_df) + 2.0),
  dpi = 300,
  bg = "white"
)

cat("GO BP terms tested:", sum(result$tested), "\n")
cat("Significant total:", nrow(sig), "\n")
cat("Significant relapse-up:", sum(sig$direction == "relapse_up"), "\n")
cat("Significant relapse-down:", sum(sig$direction == "relapse_down"), "\n")
