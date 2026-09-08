options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = "") {
  hit <- match(flag, args)
  if (!is.na(hit) && length(args) >= hit + 1) args[[hit + 1]] else default
}

default_project <- file.path(
  "C:/Users/User/Desktop/Public dataset analysis/GEO dataset/2026/single cell",
  "BTZ 相關",
  "GSE262271_TALL_relapse_scRNAseq"
)
project <- normalizePath(get_arg("--project", default_project), winslash = "/", mustWork = TRUE)
infercnv_output_name <- get_arg("--infercnv_output_name", "infercnv_output")
output_prefix <- get_arg(
  "--output_prefix",
  if (identical(infercnv_output_name, "infercnv_output")) "infercnv" else infercnv_output_name
)

paths <- list(
  objects = file.path(project, "processed_objects"),
  results = file.path(project, "results"),
  figures = file.path(project, "figures"),
  logs = file.path(project, "logs"),
  infercnv = file.path(project, "infercnv_inputs"),
  infercnv_output = file.path(project, "infercnv_inputs", infercnv_output_name)
)
invisible(lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))

prefix_file <- function(file_name) {
  if (identical(output_prefix, "infercnv")) file_name else sub("^infercnv", output_prefix, file_name)
}
result_file <- function(file_name) file.path(paths$results, prefix_file(file_name))
figure_file <- function(file_name) file.path(paths$figures, prefix_file(file_name))

user_lib <- Sys.getenv("R_LIBS_USER", unset = "C:/Users/User/Documents/R/win-library/4.4")
if (dir.exists(user_lib)) .libPaths(unique(c(normalizePath(user_lib, winslash = "/"), .libPaths())))

configure_jags <- function() {
  candidates <- unique(c(
    Sys.getenv("JAGS_HOME", unset = NA_character_),
    "C:/Users/User/Documents/JAGS/JAGS-4.3.2_extracted",
    "C:/Program Files/JAGS/JAGS-4.3.2",
    "C:/Program Files (x86)/JAGS/JAGS-4.3.2"
  ))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates) & dir.exists(candidates)]
  if (!length(candidates)) return(invisible(FALSE))
  jags_home <- normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE)
  Sys.setenv(JAGS_HOME = jags_home)
  jags_paths <- file.path(jags_home, c("x64/bin", "bin", "x64/modules", "modules"))
  jags_paths <- jags_paths[dir.exists(jags_paths)]
  if (length(jags_paths)) {
    current_path <- strsplit(Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE)[[1]]
    Sys.setenv(PATH = paste(unique(c(jags_paths, current_path)), collapse = .Platform$path.sep))
  }
  invisible(TRUE)
}
configure_jags()

load_or_stop <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
  suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))
}
load_or_stop(c("infercnv", "Seurat", "SeuratObject", "Matrix", "dplyr", "tidyr", "tibble", "ggplot2"))

write_csv_file <- function(x, path) {
  utils::write.csv(as.data.frame(x), path, row.names = FALSE, fileEncoding = "UTF-8")
}

theme_pub <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(color = "black"),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "grey70"),
      legend.key = ggplot2::element_blank()
    )
}

log_file <- file.path(paths$logs, "infercnv_burden_summary.log")
log_message <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

infercnv_file <- file.path(paths$infercnv_output, "infercnv_obj.rds")
annotated_file <- file.path(paths$objects, "GSE262271_seurat_qc_harmony_annotated.rds")
if (!file.exists(infercnv_file)) stop("Missing final inferCNV object: ", infercnv_file, call. = FALSE)
if (!file.exists(annotated_file)) stop("Missing annotated Seurat object: ", annotated_file, call. = FALSE)

log_message("Reading final inferCNV object")
infer_obj <- readRDS(infercnv_file)
expr <- infer_obj@expr.data
if (!length(dim(expr))) stop("inferCNV object has no expression matrix in @expr.data", call. = FALSE)

log_message("Computing per-cell CNV burden from denoised inferCNV expression")
if (inherits(expr, "sparseMatrix")) {
  expr@x <- pmax(expr@x, .Machine$double.eps)
  cnv_abs_log2 <- abs(log2(expr))
  cnv_abs_centered <- abs(expr - 1)
  mean_abs_log2 <- Matrix::colMeans(cnv_abs_log2, na.rm = TRUE)
  mean_abs_centered <- Matrix::colMeans(cnv_abs_centered, na.rm = TRUE)
} else {
  expr_safe <- pmax(expr, .Machine$double.eps)
  mean_abs_log2 <- colMeans(abs(log2(expr_safe)), na.rm = TRUE)
  mean_abs_centered <- colMeans(abs(expr_safe - 1), na.rm = TRUE)
}

log_message("Reading annotated Seurat metadata")
obj <- readRDS(annotated_file)
meta <- tibble::rownames_to_column(obj@meta.data, "cell")
meta_keep <- intersect(
  c("cell", "library_id", "sample_id", "hto_tag", "patient_pair", "condition", "condition_simple",
    "cluster_res0.6", "major_annotation", "broad_annotation", "infercnv_group",
    "nCount_RNA", "nFeature_RNA", "percent.mt"),
  colnames(meta)
)
meta <- meta[, meta_keep, drop = FALSE]

cell_scores <- tibble::tibble(
  cell = colnames(expr),
  infercnv_mean_abs_log2 = as.numeric(mean_abs_log2),
  infercnv_mean_abs_centered = as.numeric(mean_abs_centered)
) |>
  dplyr::left_join(meta, by = "cell") |>
  dplyr::mutate(
    infercnv_class = dplyr::case_when(
      grepl("^obs_", infercnv_group) ~ "observation_candidate_T_lineage",
      grepl("^ref_", infercnv_group) ~ "reference_non_malignant_like",
      TRUE ~ "other"
    )
  )

write_csv_file(cell_scores, result_file("infercnv_cell_cnv_burden.csv"))

summarise_burden <- function(data, group_cols) {
  data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      n_cells = dplyr::n(),
      mean_abs_log2 = mean(infercnv_mean_abs_log2, na.rm = TRUE),
      median_abs_log2 = median(infercnv_mean_abs_log2, na.rm = TRUE),
      q25_abs_log2 = unname(stats::quantile(infercnv_mean_abs_log2, 0.25, na.rm = TRUE)),
      q75_abs_log2 = unname(stats::quantile(infercnv_mean_abs_log2, 0.75, na.rm = TRUE)),
      mean_abs_centered = mean(infercnv_mean_abs_centered, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(mean_abs_log2))
}

group_summary <- summarise_burden(cell_scores, c("infercnv_class", "infercnv_group", "broad_annotation"))
sample_summary <- summarise_burden(
  cell_scores,
  c("patient_pair", "sample_id", "condition_simple", "condition", "infercnv_class", "broad_annotation")
)
cluster_summary <- summarise_burden(
  cell_scores,
  c("cluster_res0.6", "major_annotation", "broad_annotation", "infercnv_class")
)

write_csv_file(group_summary, result_file("infercnv_cnv_burden_by_group.csv"))
write_csv_file(sample_summary, result_file("infercnv_cnv_burden_by_sample.csv"))
write_csv_file(cluster_summary, result_file("infercnv_cnv_burden_by_cluster.csv"))

candidate_sample <- sample_summary |>
  dplyr::filter(broad_annotation == "T_lineage_candidate_malignant") |>
  dplyr::select(patient_pair, sample_id, condition_simple, condition, n_cells, mean_abs_log2, median_abs_log2, mean_abs_centered)
write_csv_file(candidate_sample, result_file("infercnv_candidate_malignant_like_sample_summary.csv"))

paired_candidate <- candidate_sample |>
  dplyr::select(patient_pair, condition_simple, mean_abs_log2, median_abs_log2, mean_abs_centered) |>
  tidyr::pivot_wider(
    names_from = condition_simple,
    values_from = c(mean_abs_log2, median_abs_log2, mean_abs_centered)
  ) |>
  dplyr::mutate(
    delta_relapse_minus_diagnosis_mean_abs_log2 = mean_abs_log2_relapse - mean_abs_log2_diagnosis,
    delta_relapse_minus_diagnosis_median_abs_log2 = median_abs_log2_relapse - median_abs_log2_diagnosis,
    delta_relapse_minus_diagnosis_mean_abs_centered = mean_abs_centered_relapse - mean_abs_centered_diagnosis
  )
write_csv_file(paired_candidate, result_file("infercnv_candidate_malignant_like_paired_deltas.csv"))

sample_metric_stats <- function(df, value_col) {
  wide <- df |>
    dplyr::select(patient_pair, condition_simple, value = dplyr::all_of(value_col)) |>
    tidyr::pivot_wider(names_from = condition_simple, values_from = value) |>
    dplyr::filter(!is.na(diagnosis), !is.na(relapse))
  if (nrow(wide) < 2) {
    return(tibble::tibble(metric = value_col, n_pairs = nrow(wide), mean_diagnosis = NA_real_,
                          mean_relapse = NA_real_, mean_delta = NA_real_,
                          n_pairs_relapse_higher = NA_integer_, wilcox_p = NA_real_, paired_t_p = NA_real_))
  }
  tibble::tibble(
    metric = value_col,
    n_pairs = nrow(wide),
    mean_diagnosis = mean(wide$diagnosis),
    mean_relapse = mean(wide$relapse),
    mean_delta = mean(wide$relapse - wide$diagnosis),
    n_pairs_relapse_higher = sum(wide$relapse > wide$diagnosis),
    wilcox_p = suppressWarnings(stats::wilcox.test(wide$relapse, wide$diagnosis, paired = TRUE, exact = FALSE)$p.value),
    paired_t_p = suppressWarnings(stats::t.test(wide$relapse, wide$diagnosis, paired = TRUE)$p.value)
  )
}
paired_stats <- dplyr::bind_rows(
  sample_metric_stats(candidate_sample, "mean_abs_log2"),
  sample_metric_stats(candidate_sample, "median_abs_log2"),
  sample_metric_stats(candidate_sample, "mean_abs_centered")
)
write_csv_file(paired_stats, result_file("infercnv_candidate_malignant_like_paired_statistics.csv"))

obs_ref_summary <- cell_scores |>
  dplyr::filter(infercnv_class %in% c("observation_candidate_T_lineage", "reference_non_malignant_like")) |>
  dplyr::group_by(infercnv_class) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    mean_abs_log2 = mean(infercnv_mean_abs_log2, na.rm = TRUE),
    median_abs_log2 = median(infercnv_mean_abs_log2, na.rm = TRUE),
    mean_abs_centered = mean(infercnv_mean_abs_centered, na.rm = TRUE),
    .groups = "drop"
  )
write_csv_file(obs_ref_summary, result_file("infercnv_observation_vs_reference_summary.csv"))

ordered_groups <- group_summary |>
  dplyr::mutate(label = paste0(infercnv_group, "\n", broad_annotation)) |>
  dplyr::arrange(mean_abs_log2) |>
  dplyr::pull(label)
plot_data <- cell_scores |>
  dplyr::mutate(label = paste0(infercnv_group, "\n", broad_annotation))

p_group <- ggplot2::ggplot(plot_data, ggplot2::aes(x = factor(label, levels = ordered_groups), y = infercnv_mean_abs_log2, fill = infercnv_class)) +
  ggplot2::geom_violin(scale = "width", trim = TRUE, linewidth = 0.2, alpha = 0.8) +
  ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.25, fill = "white", alpha = 0.75) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values = c(
    observation_candidate_T_lineage = "#B2433F",
    reference_non_malignant_like = "#3C6E71",
    other = "grey70"
  )) +
  ggplot2::labs(x = NULL, y = "Mean absolute log2 inferCNV deviation per cell", fill = NULL) +
  theme_pub(9)
ggplot2::ggsave(figure_file("infercnv_cnv_burden_by_group.png"), p_group, width = 8.5, height = 4.8, dpi = 300)

p_sample <- ggplot2::ggplot(candidate_sample, ggplot2::aes(x = condition_simple, y = mean_abs_log2, group = patient_pair, color = patient_pair)) +
  ggplot2::geom_line(linewidth = 0.5) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::scale_x_discrete(limits = c("diagnosis", "relapse")) +
  ggplot2::labs(x = NULL, y = "Candidate malignant-like T cells\nmean absolute log2 inferCNV deviation", color = "Pair") +
  theme_pub(10)
ggplot2::ggsave(figure_file("infercnv_candidate_malignant_like_paired_burden.png"), p_sample, width = 4.5, height = 3.6, dpi = 300)

cluster_order <- cluster_summary |>
  dplyr::arrange(mean_abs_log2) |>
  dplyr::mutate(label = paste0("C", cluster_res0.6, ": ", major_annotation)) |>
  dplyr::pull(label)
cluster_plot_data <- cell_scores |>
  dplyr::mutate(label = paste0("C", cluster_res0.6, ": ", major_annotation))
p_cluster <- ggplot2::ggplot(cluster_plot_data, ggplot2::aes(x = factor(label, levels = cluster_order), y = infercnv_mean_abs_log2, fill = broad_annotation)) +
  ggplot2::geom_boxplot(outlier.shape = NA, linewidth = 0.25) +
  ggplot2::coord_flip() +
  ggplot2::labs(x = NULL, y = "Mean absolute log2 inferCNV deviation per cell", fill = NULL) +
  theme_pub(8)
ggplot2::ggsave(figure_file("infercnv_cnv_burden_by_cluster.png"), p_cluster, width = 8.5, height = 5.8, dpi = 300)

obs_mean <- obs_ref_summary$mean_abs_log2[obs_ref_summary$infercnv_class == "observation_candidate_T_lineage"]
ref_mean <- obs_ref_summary$mean_abs_log2[obs_ref_summary$infercnv_class == "reference_non_malignant_like"]
obs_median <- obs_ref_summary$median_abs_log2[obs_ref_summary$infercnv_class == "observation_candidate_T_lineage"]
ref_median <- obs_ref_summary$median_abs_log2[obs_ref_summary$infercnv_class == "reference_non_malignant_like"]

obs_n <- sum(grepl("^obs_", cell_scores$infercnv_group), na.rm = TRUE)
ref_n <- sum(grepl("^ref_", cell_scores$infercnv_group), na.rm = TRUE)
genes_n <- nrow(expr)
cells_n <- ncol(expr)
ref_groups <- names(infer_obj@reference_grouped_cell_indices)

summary_lines <- c(
  paste0("# inferCNV summary for GSE262271 (", output_prefix, ")"),
  "",
  "## Run settings",
  "- Package: infercnv",
  "- Input: post-QC Seurat RNA counts exported from the manually annotated object",
  paste0("- Output directory: infercnv_inputs/", infercnv_output_name),
  "- Observation group: obs_T_lineage_candidate_malignant",
  paste0("- Reference groups: ", paste(ref_groups, collapse = ", ")),
  "- inferCNV run parameters: cutoff = 0.1, cluster_by_groups = TRUE, denoise = TRUE, HMM = FALSE",
  paste0("- Final denoised matrix used for burden scoring: ", format(genes_n, big.mark = ","), " genes x ", format(cells_n, big.mark = ","), " cells"),
  paste0("- Observation cells: ", format(obs_n, big.mark = ","), "; reference-like cells: ", format(ref_n, big.mark = ",")),
  "",
  "## CNV burden result",
  paste0("- Candidate T-lineage observation cells mean/median burden: ",
         signif(obs_mean, 4), " / ", signif(obs_median, 4),
         " mean absolute log2 inferCNV deviation."),
  paste0("- Reference-like cells mean/median burden: ",
         signif(ref_mean, 4), " / ", signif(ref_median, 4),
         " mean absolute log2 inferCNV deviation."),
  "- Use this as supportive evidence for candidate malignant-like T-lineage cells together with marker annotation, UMAP structure, and broad segment-level deviations in the heatmap; inspect per-reference-group burden before making a strong malignant call.",
  "- The diagnosis-relapse comparison remains paired at the patient level (N = 3 pairs); treat inferCNV burden and pathway statistics as exploratory, not definitive.",
  "",
  "## Output files",
  paste0("- infercnv_inputs/", infercnv_output_name, "/infercnv.png"),
  paste0("- figures/", prefix_file("infercnv_cnv_burden_by_group.png")),
  paste0("- figures/", prefix_file("infercnv_candidate_malignant_like_paired_burden.png")),
  paste0("- figures/", prefix_file("infercnv_cnv_burden_by_cluster.png")),
  paste0("- results/", prefix_file("infercnv_cell_cnv_burden.csv")),
  paste0("- results/", prefix_file("infercnv_cnv_burden_by_group.csv")),
  paste0("- results/", prefix_file("infercnv_cnv_burden_by_sample.csv")),
  paste0("- results/", prefix_file("infercnv_candidate_malignant_like_paired_statistics.csv"))
)
writeLines(summary_lines, result_file("infercnv_analysis_summary.md"), useBytes = TRUE)

capture.output(sessionInfo(), file = file.path(paths$logs, paste0("sessionInfo_infercnv_burden_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")))
log_message("inferCNV burden summary complete")
