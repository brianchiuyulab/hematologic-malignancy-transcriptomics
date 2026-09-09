options(stringsAsFactors = FALSE)

# Publication-style visualization of the within-sample CNV-high definition:
# (A) per-sample CNV burden distributions and the 85th-percentile cutoff,
# (B) full-KEGG GSEA robustness across CNV-high fractions, and
# (C) minimum retained cells per sample across the same fractions.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = "") {
  hit <- match(flag, args)
  if (!is.na(hit) && length(args) >= hit + 1) args[[hit + 1]] else default
}

default_project <- if (file.exists(file.path("02_Tables", "04_Method_Inputs", "infercnv_nonTref_cell_cnv_burden.csv"))) {
  "."
} else {
  file.path(
    "C:/Users/User/Desktop/Public dataset analysis/GEO dataset/2026/single cell",
    "BTZ 相關", "GSE262271_TALL_relapse_scRNAseq"
  )
}
project_arg <- get_arg("--project", default_project)
project <- if (identical(project_arg, ".")) "." else normalizePath(project_arg, winslash = "/", mustWork = TRUE)

paths <- list(
  cnv = file.path(project, "02_Tables", "04_Method_Inputs", "infercnv_nonTref_cell_cnv_burden.csv"),
  gsea = file.path(project, "02_Tables", "03_Sensitivity_Checks", "publication_compatible_candidate_full_KEGG_target_summary.csv"),
  table = file.path(project, "02_Tables", "03_Sensitivity_Checks", "within_sample_CNV_fraction_cutoff_summary.csv"),
  png = file.path(project, "01_Figures", "03_Sensitivity_Checks", "FigSens10_CNV_fraction_definition_and_pathway_robustness.png"),
  pdf = file.path(project, "01_Figures", "03_Sensitivity_Checks", "FigSens10_CNV_fraction_definition_and_pathway_robustness.pdf")
)
stopifnot(file.exists(paths$cnv), file.exists(paths$gsea))
dir.create(dirname(paths$table), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(paths$png), recursive = TRUE, showWarnings = FALSE)

user_lib <- Sys.getenv("R_LIBS_USER", unset = "C:/Users/User/Documents/R/win-library/4.4")
if (dir.exists(user_lib)) .libPaths(unique(c(normalizePath(user_lib, winslash = "/"), .libPaths())))

pkgs <- c("dplyr", "tidyr", "tibble", "ggplot2", "patchwork")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))

theme_pub <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "Arial", color = "black"),
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(color = "black"),
      strip.background = ggplot2::element_rect(fill = "grey94", color = "grey70", linewidth = 0.35),
      strip.text = ggplot2::element_text(face = "bold", size = 8.5),
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle = ggplot2::element_text(size = 8.5, color = "grey25"),
      plot.title.position = "plot",
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank()
    )
}

cnv_all <- utils::read.csv(paths$cnv, check.names = FALSE)
reference_q95 <- cnv_all |>
  dplyr::filter(grepl("^ref_", infercnv_group), is.finite(infercnv_mean_abs_log2)) |>
  dplyr::summarise(q95 = unname(stats::quantile(infercnv_mean_abs_log2, 0.95, type = 7))) |>
  dplyr::pull(q95)

cnv <- cnv_all |>
  dplyr::filter(
    broad_annotation == "T_lineage_candidate_malignant",
    is.finite(infercnv_mean_abs_log2)
  ) |>
  dplyr::mutate(
    condition_label = dplyr::recode(condition_simple, diagnosis = "Diagnosis", relapse = "Relapse"),
    sample_label = paste0(patient_pair, "  ", condition_label, "  (", sample_id, ")")
  )
candidate_fraction_above_reference_q95 <- mean(cnv$infercnv_mean_abs_log2 >= reference_q95)

sample_order <- cnv |>
  dplyr::distinct(patient_pair, condition_simple, sample_label) |>
  dplyr::arrange(patient_pair, factor(condition_simple, levels = c("diagnosis", "relapse"))) |>
  dplyr::pull(sample_label)
cnv$sample_label <- factor(cnv$sample_label, levels = sample_order)

thresholds <- cnv |>
  dplyr::group_by(sample_label, sample_id, patient_pair, condition_label) |>
  dplyr::summarise(
    candidate_cells_n = dplyr::n(),
    q85 = unname(stats::quantile(infercnv_mean_abs_log2, 0.85, type = 7)),
    selected_top15_n = sum(infercnv_mean_abs_log2 >= q85),
    .groups = "drop"
  )

x_max <- unname(stats::quantile(cnv$infercnv_mean_abs_log2, 0.999, type = 7))
density_data <- cnv |>
  dplyr::group_split(sample_label) |>
  lapply(function(z) {
    den <- stats::density(z$infercnv_mean_abs_log2, from = 0, to = x_max, n = 512, adjust = 1.0)
    tibble::tibble(sample_label = z$sample_label[[1]], x = den$x, density = den$y)
  }) |>
  dplyr::bind_rows() |>
  dplyr::left_join(thresholds |> dplyr::select(sample_label, q85), by = "sample_label")

p_a <- ggplot2::ggplot(density_data, ggplot2::aes(x = x, y = density)) +
  ggplot2::geom_area(fill = "grey88", color = NA) +
  ggplot2::geom_area(
    data = dplyr::filter(density_data, x >= q85),
    fill = "#B2182B", alpha = 0.72, color = NA
  ) +
  ggplot2::geom_line(color = "grey20", linewidth = 0.45) +
  ggplot2::geom_vline(
    data = thresholds,
    ggplot2::aes(xintercept = q85),
    color = "#B2182B", linetype = "22", linewidth = 0.55
  ) +
  ggplot2::geom_vline(
    xintercept = reference_q95,
    color = "#0072B2", linetype = "13", linewidth = 0.55
  ) +
  ggplot2::geom_text(
    data = thresholds,
    ggplot2::aes(x = q85, y = Inf, label = sprintf("q85 = %.3f", q85)),
    inherit.aes = FALSE, hjust = -0.08, vjust = 1.25, size = 2.7, color = "#8E1425"
  ) +
  ggplot2::facet_wrap(~sample_label, nrow = 2) +
  ggplot2::coord_cartesian(xlim = c(0, x_max), clip = "on") +
  ggplot2::labs(
    title = "A  Within-sample definition of the CNV-high core",
    subtitle = sprintf(
      "Red: upper 15%% per sample; red dashed: q85; blue dotted: pooled non-T reference q95 = %.3f (%.1f%% of candidate cells above)",
      reference_q95, 100 * candidate_fraction_above_reference_q95
    ),
    x = "Mean absolute inferCNV log2 deviation",
    y = "Density"
  ) +
  theme_pub(9.5) +
  ggplot2::theme(
    panel.spacing = grid::unit(0.8, "lines"),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank()
  )

fraction_ids <- c(A11 = 10, A05 = 15, A04 = 20, A02 = 25, A12 = 30, A13 = 40, A14 = 50)
gsea <- utils::read.csv(paths$gsea, check.names = FALSE) |>
  dplyr::filter(candidate_id %in% names(fraction_ids)) |>
  dplyr::mutate(
    cnv_high_percent = unname(fraction_ids[candidate_id]),
    pathway_label = dplyr::recode(pathway_name, "Mitophagy - animal" = "Mitophagy", "Proteasome" = "Proteasome"),
    FDR_tier = dplyr::case_when(
      padj < 0.05 ~ "FDR < 0.05",
      padj < 0.10 ~ "0.05 <= FDR < 0.10",
      TRUE ~ "FDR >= 0.10"
    ),
    FDR_tier = factor(FDR_tier, levels = c("FDR < 0.05", "0.05 <= FDR < 0.10", "FDR >= 0.10"))
  ) |>
  dplyr::arrange(pathway_label, cnv_high_percent)

pathway_colors <- c(Proteasome = "#0072B2", Mitophagy = "#D55E00")
p_b <- ggplot2::ggplot(
  gsea,
  ggplot2::aes(x = cnv_high_percent, y = NES, color = pathway_label, group = pathway_label)
) +
  ggplot2::annotate("rect", xmin = 13.5, xmax = 16.5, ymin = -Inf, ymax = Inf, fill = "grey70", alpha = 0.18) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_point(
    ggplot2::aes(shape = FDR_tier, fill = pathway_label),
    size = 3.0, stroke = 0.85
  ) +
  ggplot2::geom_text(
    data = dplyr::filter(gsea, cnv_high_percent == 15),
    ggplot2::aes(label = paste0("FDR ", formatC(padj, format = "f", digits = 3))),
    hjust = -0.12,
    vjust = c(-0.9, 1.5),
    size = 2.7,
    show.legend = FALSE
  ) +
  ggplot2::scale_color_manual(values = pathway_colors) +
  ggplot2::scale_fill_manual(values = pathway_colors) +
  ggplot2::scale_shape_manual(values = c("FDR < 0.05" = 16, "0.05 <= FDR < 0.10" = 17, "FDR >= 0.10" = 1)) +
  ggplot2::scale_x_continuous(breaks = sort(unique(gsea$cnv_high_percent))) +
  ggplot2::coord_cartesian(ylim = c(1.35, 2.02)) +
  ggplot2::labs(
    title = "B  Pathway robustness across CNV-high fractions",
    subtitle = "Paired pseudobulk GSEA; BH FDR across the complete KEGG collection",
    x = "CNV-high fraction within each sample (%)",
    y = "Normalized enrichment score",
    color = NULL,
    fill = NULL,
    shape = NULL
  ) +
  theme_pub(9.5) +
  ggplot2::theme(legend.position = "bottom", legend.box = "vertical")

fractions <- sort(unique(gsea$cnv_high_percent)) / 100
retention <- lapply(fractions, function(frac) {
  per_sample <- cnv |>
    dplyr::group_by(sample_id) |>
    dplyr::summarise(
      selected_n = sum(infercnv_mean_abs_log2 >= stats::quantile(
        infercnv_mean_abs_log2, 1 - frac, type = 7
      )),
      .groups = "drop"
    )
  tibble::tibble(
    cnv_high_percent = 100 * frac,
    total_cells = sum(per_sample$selected_n),
    minimum_cells_per_sample = min(per_sample$selected_n),
    maximum_cells_per_sample = max(per_sample$selected_n)
  )
}) |>
  dplyr::bind_rows()

p_c <- ggplot2::ggplot(retention, ggplot2::aes(x = cnv_high_percent, y = minimum_cells_per_sample)) +
  ggplot2::annotate("rect", xmin = 13.5, xmax = 16.5, ymin = -Inf, ymax = Inf, fill = "grey70", alpha = 0.18) +
  ggplot2::geom_hline(yintercept = 100, linetype = "22", color = "grey45", linewidth = 0.5) +
  ggplot2::geom_line(color = "grey25", linewidth = 0.7) +
  ggplot2::geom_point(shape = 21, fill = "white", color = "grey20", size = 3, stroke = 0.8) +
  ggplot2::geom_point(
    data = dplyr::filter(retention, cnv_high_percent == 15),
    shape = 21, fill = "#B2182B", color = "#B2182B", size = 3.4
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = minimum_cells_per_sample),
    vjust = -0.75, size = 2.8
  ) +
  ggplot2::annotate("text", x = 47, y = 100, label = "100 cells/sample", hjust = 1, vjust = -0.5, size = 2.7, color = "grey35") +
  ggplot2::scale_x_continuous(breaks = 100 * fractions) +
  ggplot2::expand_limits(y = max(retention$minimum_cells_per_sample) * 1.09) +
  ggplot2::labs(
    title = "C  Cell retention",
    subtitle = "Minimum retained cells among the six samples",
    x = "CNV-high fraction within each sample (%)",
    y = "Minimum cells per sample"
  ) +
  theme_pub(9.5)

cutoff_summary <- retention |>
  dplyr::left_join(
    gsea |>
      dplyr::select(cnv_high_percent, pathway_label, NES, nominal_P = pval, global_KEGG_FDR = padj) |>
      tidyr::pivot_wider(
        names_from = pathway_label,
        values_from = c(NES, nominal_P, global_KEGG_FDR),
        names_glue = "{pathway_label}_{.value}"
      ),
    by = "cnv_high_percent"
  ) |>
  dplyr::mutate(
    pooled_nonT_reference_q95 = reference_q95,
    candidate_cells_above_reference_q95_percent = 100 * candidate_fraction_above_reference_q95
  ) |>
  dplyr::arrange(cnv_high_percent)
utils::write.csv(cutoff_summary, paths$table, row.names = FALSE, fileEncoding = "UTF-8", na = "")

figure <- p_a / (p_b | p_c) +
  patchwork::plot_layout(heights = c(1.35, 1), widths = c(1.45, 1)) +
  patchwork::plot_annotation(
    title = "Defining a high-confidence CNV-high T-cell population",
    subtitle = "Source-paper QC retained; thresholds are calculated independently within each sample"
  )

png_tmp <- file.path(tempdir(), "FigSens10_CNV_fraction_definition_and_pathway_robustness.png")
ggplot2::ggsave(png_tmp, figure, width = 12.2, height = 8.8, units = "in", dpi = 320, bg = "white")
if (!file.copy(png_tmp, paths$png, overwrite = TRUE)) stop("Could not copy PNG to final path")
ggplot2::ggsave(paths$pdf, figure, width = 12.2, height = 8.8, units = "in", device = grDevices::cairo_pdf, bg = "white")
cat("Wrote:\n", paths$png, "\n", paths$pdf, "\n", paths$table, "\n", sep = "")
