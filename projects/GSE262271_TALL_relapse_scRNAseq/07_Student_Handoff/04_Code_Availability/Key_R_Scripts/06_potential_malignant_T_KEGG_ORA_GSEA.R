options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = "") {
  hit <- match(flag, args)
  if (!is.na(hit) && length(args) >= hit + 1) args[[hit + 1]] else default
}

output_root <- normalizePath(
  get_arg("--output-root", "outputs/GSE262271_TALL_relapse_scRNAseq_final"),
  winslash = "/",
  mustWork = TRUE
)
user_lib <- Sys.getenv("R_LIBS_USER", unset = "C:/Users/User/Documents/R/win-library/4.4")
if (dir.exists(user_lib)) .libPaths(unique(c(normalizePath(user_lib, winslash = "/"), .libPaths())))
required <- c("dplyr", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "), call. = FALSE)

csv_dir <- file.path(output_root, "csv", "potential_malignant_T_cells")
png_dir <- file.path(output_root, "png", "potential_malignant_T_cells")
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)

deg <- utils::read.csv(
  file.path(csv_dir, "pseudobulk_DEG_relapse_vs_diagnosis_all_genes.csv"),
  check.names = FALSE
)
kegg <- utils::read.csv(
  file.path(csv_dir, "KEGG_GSEA_all_pathways_relapse_vs_diagnosis.csv"),
  check.names = FALSE
)

deg <- deg[is.finite(deg$logFC) & is.finite(deg$PValue) & is.finite(deg$FDR) & !is.na(deg$gene), , drop = FALSE]
background <- unique(deg$gene)
background_n <- length(background)

split_genes <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  unique(unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE))
}

pathway_rows <- kegg[!duplicated(kegg$pathway_id) & kegg$tested %in% c(TRUE, "TRUE"), , drop = FALSE]
pathway_sets <- lapply(pathway_rows$linked_genes, function(x) intersect(split_genes(x), background))
names(pathway_sets) <- pathway_rows$pathway_id

run_ora <- function(selected_genes, direction, threshold_label) {
  selected_genes <- intersect(unique(selected_genes), background)
  n <- length(selected_genes)
  rows <- lapply(seq_len(nrow(pathway_rows)), function(i) {
    genes <- pathway_sets[[pathway_rows$pathway_id[[i]]]]
    K <- length(genes)
    overlap <- intersect(selected_genes, genes)
    k <- length(overlap)
    p <- if (n > 0 && K > 0 && k > 0) stats::phyper(k - 1, K, background_n - K, n, lower.tail = FALSE) else 1
    expected <- if (background_n > 0) n * K / background_n else NA_real_
    odds_ratio <- if (n > 0 && K > 0) {
      ((k + 0.5) * (background_n - K - n + k + 0.5)) /
        ((n - k + 0.5) * (K - k + 0.5))
    } else NA_real_
    data.frame(
      analysis_population = "potential malignant T cells (CNV-high top 25%)",
      comparison = "relapse_vs_diagnosis_paired_pseudobulk_N3",
      threshold = threshold_label,
      direction = direction,
      pathway_id = pathway_rows$pathway_id[[i]],
      pathway_name = pathway_rows$pathway_name[[i]],
      background_genes_n = background_n,
      selected_genes_n = n,
      pathway_genes_in_background_n = K,
      overlap_n = k,
      expected_overlap_n = expected,
      odds_ratio = odds_ratio,
      p_value = p,
      overlap_genes = paste(overlap, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$FDR <- stats::p.adjust(out$p_value, method = "BH")
  out$significant_FDR_0.05 <- out$FDR < 0.05
  out[order(out$FDR, -out$overlap_n, -out$odds_ratio), , drop = FALSE]
}

fc_cut <- log2(1.5)
strict_up <- deg$gene[deg$FDR < 0.05 & deg$logFC >= fc_cut]
strict_down <- deg$gene[deg$FDR < 0.05 & deg$logFC <= -fc_cut]
explore_up <- deg$gene[deg$PValue < 0.05 & deg$logFC >= fc_cut]
explore_down <- deg$gene[deg$PValue < 0.05 & deg$logFC <= -fc_cut]
sensitivity_up <- deg$gene[deg$PValue < 0.05 & deg$logFC >= 1.5]
sensitivity_down <- deg$gene[deg$PValue < 0.05 & deg$logFC <= -1.5]

tag_deg_set <- function(genes, label) {
  x <- deg[deg$gene %in% genes, , drop = FALSE]
  x$DEG_set <- rep(label, nrow(x))
  x
}
deg_list <- rbind(
  tag_deg_set(strict_up, "strict_FDR0.05_FC1.5_up"),
  tag_deg_set(strict_down, "strict_FDR0.05_FC1.5_down"),
  tag_deg_set(explore_up, "exploratory_nominalP0.05_FC1.5_up"),
  tag_deg_set(explore_down, "exploratory_nominalP0.05_FC1.5_down"),
  tag_deg_set(sensitivity_up, "sensitivity_nominalP0.05_abs_log2FC1.5_up"),
  tag_deg_set(sensitivity_down, "sensitivity_nominalP0.05_abs_log2FC1.5_down")
)
deg_list <- deg_list[, c("DEG_set", setdiff(names(deg_list), "DEG_set"))]
utils::write.csv(
  deg_list,
  file.path(csv_dir, "DEG_lists_for_KEGG_ORA.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

summary_table <- data.frame(
  DEG_set = c(
    "strict_FDR0.05_FC1.5_up", "strict_FDR0.05_FC1.5_down",
    "exploratory_nominalP0.05_FC1.5_up", "exploratory_nominalP0.05_FC1.5_down",
    "sensitivity_nominalP0.05_abs_log2FC1.5_up", "sensitivity_nominalP0.05_abs_log2FC1.5_down"
  ),
  definition = c(
    "FDR < 0.05 and log2FC >= log2(1.5)", "FDR < 0.05 and log2FC <= -log2(1.5)",
    "nominal P < 0.05 and log2FC >= log2(1.5)", "nominal P < 0.05 and log2FC <= -log2(1.5)",
    "nominal P < 0.05 and log2FC >= 1.5", "nominal P < 0.05 and log2FC <= -1.5"
  ),
  genes_n = c(length(strict_up), length(strict_down), length(explore_up), length(explore_down), length(sensitivity_up), length(sensitivity_down)),
  stringsAsFactors = FALSE
)
utils::write.csv(
  summary_table,
  file.path(csv_dir, "DEG_threshold_summary_for_KEGG_ORA.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

ora_exploratory <- rbind(
  run_ora(explore_up, "relapse_up", "nominal_P_lt_0.05_and_FC_ge_1.5"),
  run_ora(explore_down, "relapse_down", "nominal_P_lt_0.05_and_FC_ge_1.5")
)
ora_exploratory <- ora_exploratory[order(ora_exploratory$direction, ora_exploratory$FDR), , drop = FALSE]
utils::write.csv(
  ora_exploratory,
  file.path(csv_dir, "KEGG_ORA_exploratory_nominalP0.05_FC1.5.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

ora_sensitivity <- rbind(
  run_ora(sensitivity_up, "relapse_up", "nominal_P_lt_0.05_and_abs_log2FC_ge_1.5"),
  run_ora(sensitivity_down, "relapse_down", "nominal_P_lt_0.05_and_abs_log2FC_ge_1.5")
)
ora_sensitivity <- ora_sensitivity[order(ora_sensitivity$direction, ora_sensitivity$FDR), , drop = FALSE]
utils::write.csv(
  ora_sensitivity,
  file.path(csv_dir, "KEGG_ORA_sensitivity_nominalP0.05_abs_log2FC1.5.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

strict_status <- data.frame(
  analysis_population = "potential malignant T cells (CNV-high top 25%)",
  threshold = c("FDR < 0.05 and FC >= 1.5, relapse-up", "FDR < 0.05 and FC >= 1.5, relapse-down"),
  DEG_n = c(length(strict_up), length(strict_down)),
  KEGG_ORA_possible = c(length(strict_up) > 0, length(strict_down) > 0),
  conclusion = c(
    "No FDR-significant genes; strict KEGG ORA cannot be performed.",
    "No FDR-significant genes; strict KEGG ORA cannot be performed."
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(
  strict_status,
  file.path(csv_dir, "KEGG_ORA_strict_FDR0.05_FC1.5_status.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

gsea_significant <- kegg[kegg$tested %in% c(TRUE, "TRUE") & is.finite(kegg$padj) & kegg$padj < 0.05, , drop = FALSE]
gsea_significant <- gsea_significant[order(gsea_significant$padj, -abs(gsea_significant$NES)), , drop = FALSE]
utils::write.csv(
  gsea_significant,
  file.path(csv_dir, "KEGG_GSEA_significant_pathways_FDR0.05.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

theme_clean <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(color = "black"),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "top"
    )
}

ora_sig <- ora_exploratory[ora_exploratory$significant_FDR_0.05 & ora_exploratory$overlap_n >= 2, , drop = FALSE]
ora_candidates <- ora_exploratory[ora_exploratory$overlap_n >= 2, , drop = FALSE]
ora_plot <- do.call(rbind, lapply(split(ora_candidates, ora_candidates$direction), function(x) head(x[order(x$FDR, -x$overlap_n), ], 10)))
ora_plot$label <- paste0(ora_plot$pathway_name, " [", ifelse(ora_plot$direction == "relapse_up", "Up", "Down"), "]")
ora_plot$label <- factor(ora_plot$label, levels = rev(ora_plot$label[order(-log10(ora_plot$FDR))]))
ora_subtitle <- if (nrow(ora_sig)) {
  "DEGs: nominal P < 0.05 and absolute fold change >= 1.5; significant pathways shown"
} else {
  "No pathway reached FDR < 0.05; top exploratory enrichments are shown"
}
p_ora <- ggplot2::ggplot(
  ora_plot,
  ggplot2::aes(x = -log10(FDR), y = label, size = overlap_n, color = direction)
) +
  ggplot2::geom_vline(xintercept = -log10(0.05), linetype = 2, color = "grey45", linewidth = 0.5) +
  ggplot2::geom_point(alpha = 0.9) +
  ggplot2::scale_color_manual(values = c(relapse_up = "#B2433F", relapse_down = "#3C6E71")) +
  ggplot2::labs(
    title = "Exploratory KEGG ORA in potential malignant T cells",
    subtitle = ora_subtitle,
    x = expression(-log[10]("pathway FDR")),
    y = NULL,
    size = "Overlap genes",
    color = NULL
  ) +
  theme_clean(9)
ggplot2::ggsave(
  file.path(png_dir, "KEGG_ORA_exploratory_nominalP0.05_FC1.5.png"),
  p_ora,
  width = 9.4,
  height = max(5.2, 0.30 * nrow(ora_plot) + 1.8),
  dpi = 300,
  bg = "white"
)

gsea_plot <- head(gsea_significant, 20)
gsea_plot$pathway_label <- factor(gsea_plot$pathway_name, levels = rev(gsea_plot$pathway_name[order(gsea_plot$NES)]))
gsea_plot$fdr_label <- sprintf("FDR %.3g", gsea_plot$padj)
p_gsea <- ggplot2::ggplot(gsea_plot, ggplot2::aes(x = NES, y = pathway_label, fill = direction)) +
  ggplot2::geom_col(width = 0.72) +
  ggplot2::geom_text(ggplot2::aes(label = fdr_label), hjust = -0.05, size = 2.8) +
  ggplot2::scale_fill_manual(values = c(relapse_up = "#B2433F", relapse_down = "#3C6E71")) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.28))) +
  ggplot2::labs(
    title = "Significant KEGG GSEA pathways in potential malignant T cells",
    subtitle = "All genes ranked; relapse vs diagnosis; paired pseudobulk (N = 3); top 20 by FDR",
    x = "Normalized enrichment score (NES)",
    y = NULL,
    fill = NULL
  ) +
  theme_clean(9) +
  ggplot2::theme(legend.position = "none")
ggplot2::ggsave(
  file.path(png_dir, "KEGG_GSEA_significant_pathways_FDR0.05.png"),
  p_gsea,
  width = 9.4,
  height = 7.8,
  dpi = 300,
  bg = "white"
)

gsea_dot <- gsea_significant
gsea_dot$minus_log10_FDR <- -log10(gsea_dot$padj)
gsea_dot$pathway_label <- factor(
  gsea_dot$pathway_name,
  levels = rev(gsea_dot$pathway_name[order(gsea_dot$NES, gsea_dot$padj)])
)
p_gsea_dot <- ggplot2::ggplot(
  gsea_dot,
  ggplot2::aes(x = NES, y = pathway_label, size = size, color = minus_log10_FDR)
) +
  ggplot2::geom_point(alpha = 0.9) +
  ggplot2::scale_color_gradient(low = "#4C78A8", high = "#B2433F") +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.04, 0.08))) +
  ggplot2::labs(
    title = "Significant KEGG GSEA pathways in potential malignant T cells",
    subtitle = "Relapse vs diagnosis; paired pseudobulk (N = 3); all pathways shown have FDR < 0.05",
    x = "Normalized enrichment score (NES)",
    y = NULL,
    size = "Gene-set size",
    color = expression(-log[10]("FDR"))
  ) +
  theme_clean(9) +
  ggplot2::theme(
    legend.position = "right",
    panel.grid.major.x = ggplot2::element_line(color = "grey90", linewidth = 0.3)
  )
ggplot2::ggsave(
  file.path(png_dir, "KEGG_GSEA_significant_pathways_dotplot.png"),
  p_gsea_dot,
  width = 10.2,
  height = 9.4,
  dpi = 300,
  bg = "white"
)

cat("Strict FDR DEG counts: up=", length(strict_up), ", down=", length(strict_down), "\n", sep = "")
cat("Exploratory DEG counts: up=", length(explore_up), ", down=", length(explore_down), "\n", sep = "")
cat("Exploratory ORA significant rows: ", nrow(ora_sig), "\n", sep = "")
cat("GSEA significant pathways: ", nrow(gsea_significant), "\n", sep = "")
