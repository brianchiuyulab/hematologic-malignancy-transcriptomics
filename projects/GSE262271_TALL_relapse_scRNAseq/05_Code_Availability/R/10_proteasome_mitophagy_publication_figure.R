options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(args)) args[[1]] else "outputs/GSE262271_TALL_relapse_scRNAseq_final"
output_root <- normalizePath(output_root, winslash = "/", mustWork = TRUE)

user_lib <- Sys.getenv("R_LIBS_USER", unset = "C:/Users/User/Documents/R/win-library/4.4")
if (dir.exists(user_lib)) .libPaths(unique(c(normalizePath(user_lib, winslash = "/"), .libPaths())))
required <- c("dplyr", "ggplot2", "patchwork")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "), call. = FALSE)

csv_dir <- file.path(output_root, "csv", "potential_malignant_T_cells")
png_dir <- file.path(output_root, "png", "potential_malignant_T_cells")

gsea <- utils::read.csv(file.path(csv_dir, "KEGG_GSEA_all_pathways_relapse_vs_diagnosis.csv"), check.names = FALSE)
deg <- utils::read.csv(file.path(csv_dir, "pseudobulk_DEG_relapse_vs_diagnosis_all_genes.csv"), check.names = FALSE)
expr <- utils::read.csv(file.path(csv_dir, "proteasome_mitophagy_sample_pseudobulk_logCPM.csv"), check.names = FALSE)
module_scores <- utils::read.csv(file.path(csv_dir, "custom_mitochondrial_module_sample_scores.csv"), check.names = FALSE)

target_names <- c("Mitophagy - animal", "Proteasome")
target_gsea <- gsea[gsea$pathway_name %in% target_names, , drop = FALSE]
leading <- do.call(rbind, lapply(seq_len(nrow(target_gsea)), function(i) {
  data.frame(
    pathway_name = target_gsea$pathway_name[[i]],
    pathway_NES = target_gsea$NES[[i]],
    pathway_FDR = target_gsea$padj[[i]],
    gene = unique(strsplit(target_gsea$leadingEdge[[i]], ";", fixed = TRUE)[[1]]),
    stringsAsFactors = FALSE
  )
})) |>
  dplyr::left_join(
    deg |>
      dplyr::select(gene, logFC, PValue, FDR, rank_stat),
    by = "gene"
  )

selected <- leading |>
  dplyr::group_by(pathway_name) |>
  dplyr::slice_max(order_by = rank_stat, n = 12, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    pathway_label = dplyr::case_when(
      pathway_name == "Mitophagy - animal" ~ "Mitophagy  (NES 1.44; q = 0.0998†)",
      pathway_name == "Proteasome" ~ "Proteasome  (NES 1.53; q = 0.109)",
      TRUE ~ pathway_name
    )
  )

sample_heat <- expr |>
  dplyr::select(-dplyr::any_of("pathway_label")) |>
  dplyr::inner_join(
    selected |>
      dplyr::select(pathway_name, pathway_label, gene, rank_stat),
    by = c("pathway_name", "gene")
  ) |>
  dplyr::mutate(
    column = paste0(patient_pair, "-", ifelse(condition_simple == "diagnosis", "D", "R")),
    column_type = "Individual"
  )

mean_heat <- sample_heat |>
  dplyr::group_by(pathway_name, pathway_label, gene, rank_stat, condition_simple) |>
  dplyr::summarise(z_logCPM = mean(z_logCPM), .groups = "drop") |>
  dplyr::mutate(
    column = paste0("Mean-", ifelse(condition_simple == "diagnosis", "D", "R")),
    column_type = "Group mean"
  )
heat <- dplyr::bind_rows(
  sample_heat |>
    dplyr::select(pathway_name, pathway_label, gene, rank_stat, column, column_type, z_logCPM),
  mean_heat |>
    dplyr::select(pathway_name, pathway_label, gene, rank_stat, column, column_type, z_logCPM)
) |>
  dplyr::mutate(z_plot = pmax(-2.5, pmin(2.5, z_logCPM)))

column_levels <- c("P1-D", "P1-R", "P2-D", "P2-R", "P3-D", "P3-R", "Mean-D", "Mean-R")
heat$column <- factor(heat$column, levels = column_levels)
gene_levels <- selected |>
  dplyr::arrange(pathway_label, rank_stat) |>
  dplyr::pull(gene)
heat$gene <- factor(heat$gene, levels = unique(gene_levels))
heat$pathway_label <- factor(
  heat$pathway_label,
  levels = c("Mitophagy  (NES 1.44; q = 0.0998†)", "Proteasome  (NES 1.53; q = 0.109)")
)

p_heat <- ggplot2::ggplot(heat, ggplot2::aes(x = column, y = gene, fill = z_plot)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.25) +
  ggplot2::geom_vline(xintercept = c(2.5, 4.5), color = "grey70", linewidth = 0.35) +
  ggplot2::geom_vline(xintercept = 6.5, color = "grey25", linewidth = 0.7) +
  ggplot2::facet_grid(pathway_label ~ ., scales = "free_y", space = "free_y") +
  ggplot2::scale_fill_gradient2(
    low = "#3C6E9E", mid = "white", high = "#B2433F", midpoint = 0,
    limits = c(-2.5, 2.5), name = "Z-score"
  ) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(
    axis.text = ggplot2::element_text(color = "black"),
    axis.text.x = ggplot2::element_text(face = "bold"),
    axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", hjust = 0),
    strip.background = ggplot2::element_blank(),
    panel.spacing.y = grid::unit(0.35, "lines"),
    legend.position = "right",
    plot.margin = ggplot2::margin(4, 6, 2, 4)
  )

score <- module_scores |>
  dplyr::filter(gene_set %in% c("mitophagy_core", "proteasome_core")) |>
  dplyr::mutate(
    pathway = dplyr::recode(gene_set, mitophagy_core = "Mitophagy", proteasome_core = "Proteasome"),
    condition = factor(condition_simple, levels = c("diagnosis", "relapse"), labels = c("Diagnosis", "Relapse"))
  )

paired_stats <- do.call(rbind, lapply(split(score, score$pathway), function(x) {
  x <- x[order(x$patient_pair, x$condition), ]
  d <- x$pseudobulk_mean_logCPM[x$condition == "Diagnosis"]
  r <- x$pseudobulk_mean_logCPM[x$condition == "Relapse"]
  p_value <- stats::wilcox.test(r, d, paired = TRUE, exact = TRUE, correct = FALSE)$p.value
  yrange <- range(x$pseudobulk_mean_logCPM)
  data.frame(
    pathway = unique(x$pathway),
    x = 1.5,
    y = max(yrange) + max(diff(yrange), 0.1) * 0.22,
    label = paste0("ns\npaired P = ", formatC(p_value, digits = 2, format = "f")),
    p_value = p_value,
    stringsAsFactors = FALSE
  )
}))

p_score <- ggplot2::ggplot(
  score,
  ggplot2::aes(x = condition, y = pseudobulk_mean_logCPM, group = patient_pair)
) +
  ggplot2::geom_line(color = "grey55", linewidth = 0.55) +
  ggplot2::geom_point(ggplot2::aes(fill = condition), shape = 21, size = 2.7, color = "black", stroke = 0.4) +
  ggplot2::geom_text(
    data = paired_stats,
    ggplot2::aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 3.0,
    lineheight = 0.9
  ) +
  ggplot2::facet_wrap(~ pathway, scales = "free_y", nrow = 1) +
  ggplot2::scale_fill_manual(values = c(Diagnosis = "white", Relapse = "#B2433F")) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.28))) +
  ggplot2::labs(x = NULL, y = "Mean pseudobulk logCPM") +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(
    axis.text = ggplot2::element_text(color = "black"),
    strip.text = ggplot2::element_text(face = "bold"),
    strip.background = ggplot2::element_blank(),
    legend.position = "none",
    plot.margin = ggplot2::margin(3, 12, 4, 8)
  )

combined <- (p_heat / p_score) +
  patchwork::plot_layout(heights = c(3.1, 1.15)) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = "Proteasome and mitophagy programs in potential malignant T cells",
    subtitle = "Paired diagnosis–relapse samples (N = 3)",
    caption = "Heatmap: row-scaled pseudobulk logCPM; group means are shown for visualization only. Statistics use paired patient-level samples. †q < 0.10; ns, not significant.",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 10),
      plot.caption = ggplot2::element_text(size = 8, hjust = 0)
    )
  )

ggplot2::ggsave(
  file.path(png_dir, "proteasome_mitophagy_publication_heatmap_and_paired_scores.png"),
  combined,
  width = 8.2,
  height = 10.0,
  dpi = 400,
  bg = "white"
)

utils::write.csv(
  paired_stats,
  file.path(csv_dir, "proteasome_mitophagy_paired_pathway_statistics_for_figure.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("Saved publication figure.\n")
print(paired_stats[, c("pathway", "p_value")])
