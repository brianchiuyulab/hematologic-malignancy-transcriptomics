options(stringsAsFactors = FALSE)

# Reader-facing outputs for a fixed, defensible secondary specification:
# paper-reported QC + CNV-high top 25% within each sample. The GSEA itself was
# computed by script 17 using paired patient-level pseudobulk edgeR ranks.

args <- commandArgs(trailingOnly = TRUE)
project_arg <- if (length(args)) args[[1]] else "."
project <- if (identical(project_arg, ".")) "." else normalizePath(project_arg, winslash = "/", mustWork = TRUE)
input_file <- file.path(project, "02_Tables", "03_Sensitivity_Checks", "robustness_KEGG_GSEA_all_cases.csv")
table_dir <- file.path(project, "02_Tables", "03_Sensitivity_Checks")
figure_dir <- if (length(args) >= 2 && nzchar(args[[2]])) args[[2]] else file.path(project, "01_Figures", "03_Sensitivity_Checks")
stopifnot(file.exists(input_file))
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

user_lib <- Sys.getenv("R_LIBS_USER", unset = "C:/Users/User/Documents/R/win-library/4.4")
if (dir.exists(user_lib)) .libPaths(unique(c(normalizePath(user_lib, winslash = "/"), .libPaths())))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

write_csv_file <- function(x, path) {
  utils::write.csv(as.data.frame(x), path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
}

all_cases <- utils::read.csv(input_file, check.names = FALSE)
fixed <- all_cases |>
  dplyr::filter(case_id == "cnv_sample_top25") |>
  dplyr::arrange(padj, dplyr::desc(abs(NES)))
stopifnot(nrow(fixed) > 300, all(is.finite(fixed$NES)), all(is.finite(fixed$padj)))

fixed_fdr05 <- fixed |> dplyr::filter(padj < .05)
fixed_fdr10 <- fixed |> dplyr::filter(padj < .10)

definition <- data.frame(
  analysis_role = "Fixed secondary analysis",
  cell_population = "Candidate malignant T-lineage cells; CNV-high top 25% within each sample",
  HTO_filter = "Singlets only",
  mitochondrial_RNA_filter = "percent.mt < 5% (source-paper QC)",
  lower_feature_filter = "Sample-specific source-paper thresholds: M104/M127 >=200; M143/M148 >=300; M187/M187r >=90",
  additional_feature_floor = "None",
  upper_feature_cap = "None",
  patients_n = unique(fixed$n_patients),
  cells_n = unique(fixed$n_cells_total),
  minimum_cells_in_one_sample = unique(fixed$min_cells_per_sample),
  pseudobulk_design = "~ patient_pair + condition_simple",
  GSEA_rank = "sign(logFC) * sqrt(F)",
  pathway_collection = "Complete KEGG collection used in the primary analysis",
  multiple_testing = "BH FDR across all tested KEGG pathways",
  stringsAsFactors = FALSE
)

mechanism_map <- data.frame(
  pathway_id = c(
    "hsa03010", "hsa03008", "hsa03040",
    "hsa00190", "hsa05208", "hsa03050", "hsa04137",
    "hsa04658", "hsa04659", "hsa04630", "hsa04660", "hsa04060", "hsa05235",
    "hsa04010", "hsa04064", "hsa04066", "hsa04210", "hsa04217", "hsa04218", "hsa05202",
    "hsa00511", "hsa04742", "hsa00562", "hsa04070"
  ),
  display_name = c(
    "Ribosome", "Ribosome biogenesis", "Spliceosome",
    "Oxidative phosphorylation", "Reactive oxygen species program", "Proteasome", "Mitophagy",
    "Th1 and Th2 differentiation", "Th17 differentiation", "JAK-STAT signaling", "T-cell receptor signaling",
    "Cytokine-receptor interaction", "PD-L1 and PD-1 checkpoint",
    "MAPK signaling", "NF-kappa B signaling", "HIF-1 signaling", "Apoptosis", "Necroptosis",
    "Cellular senescence", "Transcriptional misregulation in cancer",
    "Other glycan degradation", "Taste transduction", "Inositol phosphate metabolism",
    "Phosphatidylinositol signaling"
  ),
  program = c(
    rep("Translation and RNA processing", 3),
    rep("Mitochondria and proteostasis", 4),
    rep("T-cell and cytokine signaling", 6),
    rep("Stress, survival and cell fate", 7),
    rep("Relapse-down pathways", 4)
  ),
  stringsAsFactors = FALSE
)

curated <- fixed_fdr10 |>
  dplyr::inner_join(mechanism_map, by = "pathway_id") |>
  dplyr::mutate(
    program = factor(program, levels = c(
      "Translation and RNA processing", "Mitochondria and proteostasis",
      "T-cell and cytokine signaling", "Stress, survival and cell fate",
      "Relapse-down pathways"
    )),
    FDR_tier = ifelse(padj < .05, "FDR < 0.05", "FDR 0.05 to <0.10"),
    FDR_tier = factor(FDR_tier, levels = c("FDR < 0.05", "FDR 0.05 to <0.10")),
    FDR_evidence = pmin(-log10(padj), 6),
    FDR_label = ifelse(padj < .001, formatC(padj, format = "e", digits = 1), sprintf("%.3f", padj))
  ) |>
  dplyr::arrange(program, NES) |>
  dplyr::mutate(display_name = factor(display_name, levels = unique(display_name)))

write_csv_file(definition, file.path(table_dir, "fixed_secondary_within_sample_top25_definition.csv"))
write_csv_file(fixed_fdr05, file.path(table_dir, "fixed_secondary_within_sample_top25_KEGG_GSEA_FDR0.05.csv"))
write_csv_file(fixed_fdr10, file.path(table_dir, "fixed_secondary_within_sample_top25_KEGG_GSEA_FDR0.10.csv"))
write_csv_file(curated, file.path(table_dir, "fixed_secondary_within_sample_top25_KEGG_curated_mechanistic_pathways.csv"))

palette <- c(
  "Translation and RNA processing" = "#0072B2",
  "Mitochondria and proteostasis" = "#D55E00",
  "T-cell and cytokine signaling" = "#009E73",
  "Stress, survival and cell fate" = "#CC79A7",
  "Relapse-down pathways" = "#7A7A7A"
)

p <- ggplot2::ggplot(curated, ggplot2::aes(y = display_name, x = NES)) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = NES, yend = display_name),
    color = "grey78", linewidth = .55
  ) +
  ggplot2::geom_point(
    ggplot2::aes(fill = program, shape = FDR_tier, size = FDR_evidence),
    color = "black", stroke = .45
  ) +
  ggplot2::geom_text(
    ggplot2::aes(x = 2.70, label = paste0("FDR ", FDR_label)),
    hjust = 0, size = 3.0, color = "grey20"
  ) +
  ggplot2::facet_grid(program ~ ., scales = "free_y", space = "free_y", switch = "y") +
  ggplot2::scale_fill_manual(values = palette, guide = "none") +
  ggplot2::scale_shape_manual(values = c("FDR < 0.05" = 21, "FDR 0.05 to <0.10" = 24)) +
  ggplot2::scale_size_continuous(
    range = c(3.1, 7.0), limits = c(1, 6), breaks = c(2, 4, 6),
    guide = "none"
  ) +
  ggplot2::coord_cartesian(xlim = c(-2.15, 3.22), clip = "off") +
  ggplot2::labs(
    title = "Relapse-associated KEGG programs in potential malignant T cells",
    subtitle = "CNV-high top 25% within each sample; source-paper QC; paired patient-level pseudobulk GSEA",
    x = "Normalized enrichment score (relapse vs diagnosis)", y = NULL,
    shape = "Global KEGG FDR"
  ) +
  ggplot2::theme_classic(base_size = 10.5) +
  ggplot2::theme(
    text = ggplot2::element_text(family = "Arial", color = "black"),
    axis.text = ggplot2::element_text(color = "black"),
    axis.text.y = ggplot2::element_text(size = 9.2),
    axis.line.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    strip.placement = "outside",
    strip.background = ggplot2::element_blank(),
    strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold", hjust = 1, size = 9.2),
    panel.spacing.y = grid::unit(0.55, "lines"),
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.title.position = "plot",
    plot.margin = ggplot2::margin(8, 34, 8, 8)
  )

ggplot2::ggsave(
  file.path(figure_dir, "FigSens09_within_sample_top25_KEGG_mechanistic_programs.png"),
  p, width = 10.8, height = 10.8, dpi = 360, bg = "white"
)

cat(
  "Fixed secondary KEGG: tested", nrow(fixed),
  "pathways; FDR<0.05", nrow(fixed_fdr05),
  "; FDR<0.10", nrow(fixed_fdr10),
  "; curated figure", nrow(curated), "pathways.\n"
)
