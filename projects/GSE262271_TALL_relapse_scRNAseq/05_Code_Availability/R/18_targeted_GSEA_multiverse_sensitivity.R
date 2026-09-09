options(stringsAsFactors = FALSE)

# Large targeted sensitivity grid for KEGG Proteasome and Mitophagy - animal.
# This script reports every predefined specification and separately exports all
# nominal-P < 0.05 rows. It does not replace the frozen primary KEGG analysis.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = "") {
  hit <- match(flag, args)
  if (!is.na(hit) && length(args) >= hit + 1) args[[hit + 1]] else default
}

default_project <- if (file.exists(file.path("04_R_Objects", "GSE262271_seurat_qc_harmony_annotated.rds"))) "." else file.path(
  "C:/Users/User/Desktop/Public dataset analysis/GEO dataset/2026/single cell",
  "BTZ 相關", "GSE262271_TALL_relapse_scRNAseq"
)
project_arg <- get_arg("--project", default_project)
project <- if (identical(project_arg, ".")) "." else normalizePath(project_arg, winslash = "/", mustWork = TRUE)
nperm <- as.integer(get_arg("--nperm", "20000"))
if (!is.finite(nperm) || nperm < 1000) nperm <- 20000L

paths <- list(
  object = file.path(project, "04_R_Objects", "GSE262271_seurat_qc_harmony_annotated.rds"),
  cnv = file.path(project, "02_Tables", "04_Method_Inputs", "infercnv_nonTref_cell_cnv_burden.csv"),
  kegg = file.path(project, "02_Tables", "02_Full_Results", "KEGG_GSEA_all_pathways_relapse_vs_diagnosis.csv"),
  tables = file.path(project, "02_Tables", "03_Sensitivity_Checks"),
  figures = file.path(project, "01_Figures", "03_Sensitivity_Checks"),
  environment = file.path(project, "05_Code_Availability", "environment")
)
stopifnot(file.exists(paths$object), file.exists(paths$cnv), file.exists(paths$kegg))
dir.create(paths$tables, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$figures, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$environment, recursive = TRUE, showWarnings = FALSE)

user_lib <- Sys.getenv("R_LIBS_USER", unset = "C:/Users/User/Documents/R/win-library/4.4")
if (dir.exists(user_lib)) .libPaths(unique(c(normalizePath(user_lib, winslash = "/"), .libPaths())))
pkgs <- c("Seurat", "SeuratObject", "Matrix", "edgeR", "fgsea", "dplyr", "tidyr", "tibble", "ggplot2")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))

write_csv_file <- function(x, path) {
  utils::write.csv(as.data.frame(x), path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
}

theme_pub <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "Arial", color = "black"),
      axis.text = ggplot2::element_text(color = "black"),
      strip.background = ggplot2::element_rect(fill = "grey94", color = "grey75", linewidth = 0.35),
      strip.text = ggplot2::element_text(face = "bold", color = "black"),
      legend.key = ggplot2::element_blank(),
      plot.title.position = "plot"
    )
}

cat("Loading post-QC object...\n")
obj <- readRDS(paths$object)
if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) obj <- SeuratObject::JoinLayers(obj, assay = "RNA")
counts <- Seurat::GetAssayData(obj, assay = "RNA", layer = "counts")
meta <- tibble::rownames_to_column(obj@meta.data, "cell")
cnv <- utils::read.csv(paths$cnv, check.names = FALSE)[, c("cell", "infercnv_mean_abs_log2")]
meta <- meta |>
  dplyr::select(-dplyr::any_of("infercnv_mean_abs_log2")) |>
  dplyr::left_join(cnv, by = "cell") |>
  dplyr::mutate(
    patient_pair = factor(patient_pair, levels = c("P1", "P2", "P3")),
    condition_simple = factor(condition_simple, levels = c("diagnosis", "relapse"))
  )
candidate <- meta |>
  dplyr::filter(broad_annotation == "T_lineage_candidate_malignant", is.finite(infercnv_mean_abs_log2))
sample_info <- candidate |>
  dplyr::distinct(sample_id, patient_pair, condition_simple) |>
  dplyr::arrange(patient_pair, condition_simple)

kegg <- utils::read.csv(paths$kegg, check.names = FALSE)
target_rows <- kegg |>
  dplyr::filter(pathway_id %in% c("hsa03050", "hsa04137"), !duplicated(pathway_id))
target_sets <- strsplit(target_rows$linked_genes, ";", fixed = TRUE)
target_sets <- lapply(target_sets, function(x) sort(unique(trimws(x[nzchar(trimws(x))]))))
names(target_sets) <- target_rows$pathway_name

sample_bounds <- candidate |>
  dplyr::group_by(sample_id) |>
  dplyr::summarise(nf_q99 = unname(stats::quantile(nFeature_RNA, .99)), .groups = "drop")
candidate <- candidate |> dplyr::left_join(sample_bounds, by = "sample_id")

population_grid <- tidyr::crossing(
  selection_scope = c("global", "within_sample"),
  top_fraction = c(.15, .20, .25, .30, .40, .50),
  mt_max = c(3.0, 3.5, 4.0, 4.5, 5.0),
  min_features = c(0L, 500L, 1000L, 1500L),
  upper_feature_cap = c("none", "sample_q99")
) |>
  dplyr::mutate(
    specification_axis = "Population + QC",
    gene_filter = "filterByExpr",
    robust_dispersion = TRUE,
    normalization = "TMM",
    rank_metric = "signed_sqrt_F",
    specification_id = sprintf("PQ%03d", dplyr::row_number())
  )

model_grid <- tidyr::crossing(
  gene_filter = c("filterByExpr", "cpm05_n2", "cpm1_n2", "cpm1_n3"),
  robust_dispersion = c(TRUE, FALSE),
  normalization = c("TMM", "TMMwsp"),
  rank_metric = c("signed_sqrt_F", "signed_neglog10_P", "logFC")
) |>
  dplyr::mutate(
    specification_axis = "Model + ranking",
    selection_scope = "global",
    top_fraction = .25,
    mt_max = 5.0,
    min_features = 0L,
    upper_feature_cap = "none",
    specification_id = sprintf("MR%03d", dplyr::row_number())
  )

specifications <- dplyr::bind_rows(population_grid, model_grid) |>
  dplyr::mutate(
    primary_equivalent = selection_scope == "global" & top_fraction == .25 & mt_max == 5 &
      min_features == 0 & upper_feature_cap == "none" & gene_filter == "filterByExpr" &
      robust_dispersion & normalization == "TMM" & rank_metric == "signed_sqrt_F"
  )

select_cells <- function(spec) {
  x <- candidate
  if (spec$selection_scope == "within_sample") {
    x <- x |>
      dplyr::group_by(sample_id) |>
      dplyr::mutate(threshold = unname(stats::quantile(infercnv_mean_abs_log2, 1 - spec$top_fraction))) |>
      dplyr::ungroup() |>
      dplyr::filter(infercnv_mean_abs_log2 >= threshold)
  } else {
    threshold <- unname(stats::quantile(candidate$infercnv_mean_abs_log2, 1 - spec$top_fraction))
    x <- x |> dplyr::filter(infercnv_mean_abs_log2 >= threshold)
  }
  x <- x |> dplyr::filter(percent.mt < spec$mt_max)
  if (spec$min_features > 0) x <- x |> dplyr::filter(nFeature_RNA >= spec$min_features)
  if (spec$upper_feature_cap == "sample_q99") x <- x |> dplyr::filter(nFeature_RNA <= nf_q99)
  x
}

make_pb <- function(x) {
  idx <- match(x$sample_id, sample_info$sample_id)
  map <- Matrix::sparseMatrix(
    i = seq_len(nrow(x)), j = idx, x = 1,
    dims = c(nrow(x), nrow(sample_info)),
    dimnames = list(x$cell, sample_info$sample_id)
  )
  as(counts[, x$cell, drop = FALSE] %*% map, "dgCMatrix")
}

keep_genes <- function(y, design, method) {
  if (method == "filterByExpr") return(edgeR::filterByExpr(y, design = design))
  z <- edgeR::cpm(y)
  if (method == "cpm05_n2") return(rowSums(z >= .5) >= 2)
  if (method == "cpm1_n2") return(rowSums(z >= 1) >= 2)
  if (method == "cpm1_n3") return(rowSums(z >= 1) >= 3)
  stop("Unknown filter: ", method)
}

fit_deg <- function(x, gene_filter, robust_dispersion, normalization) {
  pb <- make_pb(x)
  y <- edgeR::DGEList(counts = pb, samples = as.data.frame(sample_info))
  design <- stats::model.matrix(~ patient_pair + condition_simple, data = sample_info)
  keep <- keep_genes(y, design, gene_filter)
  if (sum(keep) < 100) return(NULL)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y, method = normalization)
  y <- edgeR::estimateDisp(y, design, robust = robust_dispersion)
  fit <- edgeR::glmQLFit(y, design, robust = robust_dispersion)
  coef_name <- grep("condition_simple.*relapse", colnames(design), value = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = coef_name)
  deg <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
  deg$gene <- rownames(deg)
  deg
}

make_rank <- function(deg, metric) {
  value <- switch(
    metric,
    signed_sqrt_F = sign(deg$logFC) * sqrt(pmax(deg$F, 0)),
    signed_neglog10_P = sign(deg$logFC) * -log10(pmax(deg$PValue, 1e-300)),
    logFC = deg$logFC,
    stop("Unknown rank metric: ", metric)
  )
  names(value) <- deg$gene
  value <- value[is.finite(value) & !duplicated(names(value))]
  value <- value + rank(names(value), ties.method = "first") * 1e-12
  sort(value, decreasing = TRUE)
}

run_target_gsea <- function(rank_stats, seed) {
  sets <- lapply(target_sets, intersect, y = names(rank_stats))
  set.seed(seed)
  fg <- suppressWarnings(fgsea::fgseaSimple(
    pathways = sets, stats = rank_stats, nperm = nperm,
    minSize = 5, maxSize = 500, nproc = 1
  ))
  fg <- as.data.frame(fg)
  fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  fg
}

summarize_selection <- function(x) {
  per_sample <- x |> dplyr::count(sample_id, name = "n_cells")
  n_by_sample <- sample_info |>
    dplyr::left_join(per_sample, by = "sample_id") |>
    dplyr::mutate(n_cells = ifelse(is.na(n_cells), 0L, n_cells))
  state_prop <- x |>
    dplyr::count(major_annotation) |>
    dplyr::mutate(prop = n / sum(n))
  get_prop <- function(label) {
    z <- state_prop$prop[state_prop$major_annotation == label]
    if (length(z)) z[[1]] else 0
  }
  tibble::tibble(
    n_cells_total = nrow(x),
    min_cells_per_sample = min(n_by_sample$n_cells),
    max_cells_per_sample = max(n_by_sample$n_cells),
    sample_cell_imbalance_ratio = max(n_by_sample$n_cells) / max(1, min(n_by_sample$n_cells)),
    diagnosis_cells_n = sum(n_by_sample$n_cells[n_by_sample$condition_simple == "diagnosis"]),
    relapse_cells_n = sum(n_by_sample$n_cells[n_by_sample$condition_simple == "relapse"]),
    median_percent_mt = stats::median(x$percent.mt),
    median_nFeature_RNA = stats::median(x$nFeature_RNA),
    cycling_prop = get_prop("Cycling_TALL_blast_like"),
    immature_prop = get_prop("Immature_HSC_like_TALL"),
    quiescent_prop = get_prop("Quiescent_CD44_high_TALL_like"),
    er_stress_prop = get_prop("ER_stress_TALL_like")
  )
}

results <- list()
registry <- list()
result_index <- 0L
registry_index <- 0L

cat("Running", nrow(population_grid), "population/QC combinations...\n")
for (i in seq_len(nrow(population_grid))) {
  if (i %% 20 == 1) cat(sprintf("Population/QC %d/%d\n", i, nrow(population_grid)))
  spec <- population_grid[i, , drop = FALSE]
  x <- select_cells(spec)
  selection_summary <- summarize_selection(x)
  valid <- selection_summary$min_cells_per_sample >= 20
  registry_index <- registry_index + 1L
  registry[[registry_index]] <- dplyr::bind_cols(spec, selection_summary) |>
    dplyr::mutate(valid = valid, invalid_reason = ifelse(valid, "", "fewer than 20 cells in at least one sample"))
  if (!valid) next
  deg <- fit_deg(x, spec$gene_filter, spec$robust_dispersion, spec$normalization)
  if (is.null(deg)) next
  ranks <- make_rank(deg, spec$rank_metric)
  fg <- run_target_gsea(ranks, 180000L + i)
  result_index <- result_index + 1L
  results[[result_index]] <- dplyr::bind_cols(spec[rep(1, nrow(fg)), , drop = FALSE], selection_summary[rep(1, nrow(fg)), , drop = FALSE], fg)
}

cat("Running", nrow(model_grid), "model/ranking combinations...\n")
primary_selection_spec <- population_grid |>
  dplyr::filter(selection_scope == "global", top_fraction == .25, mt_max == 5, min_features == 0, upper_feature_cap == "none") |>
  dplyr::slice(1)
primary_cells <- select_cells(primary_selection_spec)
primary_selection_summary <- summarize_selection(primary_cells)
model_groups <- model_grid |>
  dplyr::distinct(gene_filter, robust_dispersion, normalization)
for (g in seq_len(nrow(model_groups))) {
  grp <- model_groups[g, , drop = FALSE]
  cat(sprintf("Model group %d/%d: %s, robust=%s, %s\n", g, nrow(model_groups), grp$gene_filter, grp$robust_dispersion, grp$normalization))
  deg <- fit_deg(primary_cells, grp$gene_filter, grp$robust_dispersion, grp$normalization)
  matching <- model_grid |>
    dplyr::filter(
      gene_filter == grp$gene_filter,
      robust_dispersion == grp$robust_dispersion,
      normalization == grp$normalization
    )
  for (j in seq_len(nrow(matching))) {
    spec <- matching[j, , drop = FALSE]
    registry_index <- registry_index + 1L
    registry[[registry_index]] <- dplyr::bind_cols(spec, primary_selection_summary) |>
      dplyr::mutate(valid = !is.null(deg), invalid_reason = ifelse(is.null(deg), "fewer than 100 genes after filtering", ""))
    if (is.null(deg)) next
    ranks <- make_rank(deg, spec$rank_metric)
    fg <- run_target_gsea(ranks, 280000L + match(spec$specification_id, model_grid$specification_id))
    result_index <- result_index + 1L
    results[[result_index]] <- dplyr::bind_cols(spec[rep(1, nrow(fg)), , drop = FALSE], primary_selection_summary[rep(1, nrow(fg)), , drop = FALSE], fg)
  }
}

registry_all <- dplyr::bind_rows(registry) |>
  dplyr::mutate(
    primary_equivalent = selection_scope == "global" & top_fraction == .25 & mt_max == 5 &
      min_features == 0 & upper_feature_cap == "none" & gene_filter == "filterByExpr" &
      robust_dispersion & normalization == "TMM" & rank_metric == "signed_sqrt_F"
  )
result_all <- dplyr::bind_rows(results) |>
  dplyr::rename(pathway_name = pathway) |>
  dplyr::group_by(specification_id) |>
  dplyr::mutate(BH_within_spec_two_targets = stats::p.adjust(pval, method = "BH")) |>
  dplyr::ungroup() |>
  dplyr::group_by(pathway_name) |>
  dplyr::mutate(BH_across_specifications_same_pathway = stats::p.adjust(pval, method = "BH")) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    BH_across_all_specifications_and_targets = stats::p.adjust(pval, method = "BH"),
    nominal_P_monte_carlo_SE = sqrt(pval * (1 - pval) / nperm),
    nominal_P_borderline_0.05 = abs(pval - .05) <= 2 * nominal_P_monte_carlo_SE,
    nominal_P_lt_0.05 = pval < .05,
    two_target_BH_lt_0.05 = BH_within_spec_two_targets < .05,
    direction = ifelse(NES > 0, "relapse_up", "relapse_down"),
    nperm = nperm,
    primary_equivalent = selection_scope == "global" & top_fraction == .25 & mt_max == 5 &
      min_features == 0 & upper_feature_cap == "none" & gene_filter == "filterByExpr" &
      robust_dispersion & normalization == "TMM" & rank_metric == "signed_sqrt_F"
  ) |>
  dplyr::arrange(pathway_name, pval, specification_id)

hits <- result_all |>
  dplyr::filter(nominal_P_lt_0.05) |>
  dplyr::arrange(pathway_name, pval)

# Reader-facing extracts: the best results per pathway and one row for every
# specification in which both prespecified targets have nominal P < 0.05.
top_combinations <- result_all |>
  dplyr::group_by(pathway_name) |>
  dplyr::slice_min(order_by = pval, n = 25, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::arrange(pathway_name, pval)

specification_metadata <- result_all |>
  dplyr::distinct(specification_id, .keep_all = TRUE) |>
  dplyr::select(
    specification_id, specification_axis, selection_scope, top_fraction,
    mt_max, min_features, upper_feature_cap, gene_filter, robust_dispersion,
    normalization, rank_metric, n_cells_total, min_cells_per_sample,
    sample_cell_imbalance_ratio, diagnosis_cells_n, relapse_cells_n,
    median_percent_mt, median_nFeature_RNA, cycling_prop, immature_prop,
    quiescent_prop, er_stress_prop, nperm, primary_equivalent
  )
proteasome_wide <- result_all |>
  dplyr::filter(pathway_name == "Proteasome") |>
  dplyr::transmute(
    specification_id, proteasome_NES = NES, proteasome_nominal_P = pval,
    proteasome_two_target_BH = BH_within_spec_two_targets
  )
mitophagy_wide <- result_all |>
  dplyr::filter(pathway_name == "Mitophagy - animal") |>
  dplyr::transmute(
    specification_id, mitophagy_NES = NES, mitophagy_nominal_P = pval,
    mitophagy_two_target_BH = BH_within_spec_two_targets
  )
both_target_hits <- specification_metadata |>
  dplyr::left_join(proteasome_wide, by = "specification_id") |>
  dplyr::left_join(mitophagy_wide, by = "specification_id") |>
  dplyr::filter(proteasome_nominal_P < .05, mitophagy_nominal_P < .05) |>
  dplyr::mutate(worst_target_nominal_P = pmax(proteasome_nominal_P, mitophagy_nominal_P)) |>
  dplyr::arrange(worst_target_nominal_P, specification_id)
summary_by_pathway <- result_all |>
  dplyr::group_by(pathway_name) |>
  dplyr::summarise(
    valid_specifications_n = dplyr::n(),
    relapse_up_n = sum(NES > 0),
    nominal_P_lt_0.05_n = sum(nominal_P_lt_0.05),
    nominal_P_lt_0.05_pct = mean(nominal_P_lt_0.05) * 100,
    two_target_BH_lt_0.05_n = sum(two_target_BH_lt_0.05),
    two_target_BH_lt_0.05_pct = mean(two_target_BH_lt_0.05) * 100,
    nominal_P_min = min(pval),
    nominal_P_median = stats::median(pval),
    nominal_P_max = max(pval),
    NES_min = min(NES), NES_median = stats::median(NES), NES_max = max(NES),
    .groups = "drop"
  )

parameter_summary_one <- function(data, variable) {
  data |>
    dplyr::mutate(parameter = variable, level = as.character(.data[[variable]])) |>
    dplyr::group_by(pathway_name, specification_axis, parameter, level) |>
    dplyr::summarise(
      specifications_n = dplyr::n(),
      nominal_hit_n = sum(nominal_P_lt_0.05),
      nominal_hit_pct = mean(nominal_P_lt_0.05) * 100,
      median_nominal_P = stats::median(pval),
      minimum_nominal_P = min(pval),
      median_NES = stats::median(NES),
      median_cells_n = stats::median(n_cells_total),
      .groups = "drop"
    )
}
parameter_summary <- dplyr::bind_rows(lapply(
  c("selection_scope", "top_fraction", "mt_max", "min_features", "upper_feature_cap",
    "gene_filter", "robust_dispersion", "normalization", "rank_metric"),
  function(v) parameter_summary_one(result_all, v)
))

write_csv_file(registry_all, file.path(paths$tables, "multiverse_specification_registry.csv"))
write_csv_file(result_all, file.path(paths$tables, "multiverse_targeted_GSEA_all_results.csv"))
write_csv_file(hits, file.path(paths$tables, "multiverse_nominal_P_lt_0.05_combinations.csv"))
write_csv_file(top_combinations, file.path(paths$tables, "multiverse_top25_combinations_by_pathway.csv"))
write_csv_file(both_target_hits, file.path(paths$tables, "multiverse_both_targets_nominal_P_lt_0.05.csv"))
write_csv_file(summary_by_pathway, file.path(paths$tables, "multiverse_targeted_GSEA_summary.csv"))
write_csv_file(parameter_summary, file.path(paths$tables, "multiverse_parameter_level_summary.csv"))

# Specification curve: every valid result is shown.
curve <- result_all |>
  dplyr::group_by(pathway_name) |>
  dplyr::arrange(pval, .by_group = TRUE) |>
  dplyr::mutate(specification_rank = dplyr::row_number()) |>
  dplyr::ungroup()
curve_annot <- summary_by_pathway |>
  dplyr::mutate(
    x = Inf, y = Inf,
    label = paste0(nominal_P_lt_0.05_n, "/", valid_specifications_n, " nominal P < 0.05")
  )
p_curve <- ggplot2::ggplot(curve, ggplot2::aes(
  x = specification_rank, y = -log10(pval), color = specification_axis,
  shape = two_target_BH_lt_0.05
)) +
  ggplot2::geom_hline(yintercept = -log10(.05), linetype = 2, color = "grey45", linewidth = .5) +
  ggplot2::geom_point(alpha = .68, size = 1.7, stroke = .35) +
  ggplot2::geom_text(
    data = curve_annot, ggplot2::aes(x = x, y = y, label = label),
    inherit.aes = FALSE, hjust = 1.05, vjust = 1.4, size = 3.2
  ) +
  ggplot2::facet_wrap(~ pathway_name, ncol = 1, scales = "free_y") +
  ggplot2::scale_color_manual(values = c("Population + QC" = "#0072B2", "Model + ranking" = "#D55E00")) +
  ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "< 0.05", `FALSE` = ">= 0.05")) +
  ggplot2::labs(
    title = "Targeted GSEA multiverse specification curve",
    subtitle = paste0("All predefined valid specifications; ", format(nperm, big.mark = ","), " gene-label permutations per specification"),
    x = "Specifications ordered by nominal P within pathway",
    y = expression(-log[10]("nominal P")),
    color = "Sensitivity axis", shape = "BH across 2 targets"
  ) +
  theme_pub(10) +
  ggplot2::theme(legend.position = "bottom", panel.grid.major.x = ggplot2::element_blank())
ggplot2::ggsave(
  file.path(paths$figures, "FigSens07_targeted_GSEA_multiverse_specification_curve.png"),
  p_curve, width = 9.0, height = 7.8, dpi = 360, bg = "white"
)

# Hit-rate heatmap for the 480 population/QC combinations.
heat <- result_all |>
  dplyr::filter(specification_axis == "Population + QC") |>
  dplyr::group_by(pathway_name, selection_scope, top_fraction, mt_max) |>
  dplyr::summarise(
    hit_rate = mean(nominal_P_lt_0.05),
    valid_combinations_n = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    top_label = paste0("Top ", round(top_fraction * 100), "%"),
    top_label = factor(top_label, levels = paste0("Top ", c(15, 20, 25, 30, 40, 50), "%")),
    mt_label = factor(paste0("<", mt_max, "%"), levels = paste0("<", c(5, 4.5, 4, 3.5, 3), "%"))
  )
p_heat <- ggplot2::ggplot(heat, ggplot2::aes(top_label, mt_label, fill = hit_rate)) +
  ggplot2::geom_tile(color = "white", linewidth = .7) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.0f%%", hit_rate * 100)), size = 3.0) +
  ggplot2::facet_grid(
    pathway_name ~ selection_scope,
    labeller = ggplot2::labeller(
      selection_scope = c(global = "Global CNV threshold", within_sample = "Within-sample CNV threshold"),
      pathway_name = c(`Mitophagy - animal` = "KEGG Mitophagy", Proteasome = "KEGG Proteasome")
    )
  ) +
  ggplot2::scale_fill_gradient(low = "#F2F2F2", high = "#B2182B", limits = c(0, 1), labels = function(x) paste0(round(x * 100), "%")) +
  ggplot2::labs(
    title = "Frequency of nominal P < 0.05 across QC combinations",
    subtitle = "Each tile summarizes min-feature and upper-feature-cap combinations",
    x = "CNV-high fraction", y = "Mitochondrial RNA threshold", fill = "Hit rate"
  ) +
  theme_pub(9) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
    axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
    strip.text.y = ggplot2::element_text(angle = 0)
  )
ggplot2::ggsave(
  file.path(paths$figures, "FigSens08_multiverse_QC_hit_rate_heatmap.png"),
  p_heat, width = 9.2, height = 7.2, dpi = 360, bg = "white"
)

capture.output(sessionInfo(), file = file.path(paths$environment, "sessionInfo_18_targeted_GSEA_multiverse_sensitivity.txt"))
cat("Completed", nrow(registry_all), "specifications and", nrow(result_all), "pathway results.\n")
