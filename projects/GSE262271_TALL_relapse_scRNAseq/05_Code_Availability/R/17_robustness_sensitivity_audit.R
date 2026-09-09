options(stringsAsFactors = FALSE)
options(timeout = 300)

# Robustness audit for the paired GSE262271 relapse analysis.
# The primary analysis remains fixed: global CNV-high top 25% among candidate
# malignant T-lineage cells, original post-QC object, paired edgeR pseudobulk,
# filterByExpr, robust QL dispersion, and sign(logFC) * sqrt(F) GSEA ranking.
# Alternative settings below are one-factor-at-a-time sensitivity analyses.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = "") {
  hit <- match(flag, args)
  if (!is.na(hit) && length(args) >= hit + 1) args[[hit + 1]] else default
}

default_project <- if (file.exists(file.path("04_R_Objects", "GSE262271_seurat_qc_harmony_annotated.rds"))) {
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
  object = file.path(project, "04_R_Objects", "GSE262271_seurat_qc_harmony_annotated.rds"),
  cnv = file.path(project, "02_Tables", "04_Method_Inputs", "infercnv_nonTref_cell_cnv_burden.csv"),
  kegg = file.path(project, "02_Tables", "02_Full_Results", "KEGG_GSEA_all_pathways_relapse_vs_diagnosis.csv"),
  cross = file.path(project, "02_Tables", "01_Key_Results", "cross_dataset_KEGG_FDR0.10_same_direction_intersection.csv"),
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
      axis.title = ggplot2::element_text(color = "black"),
      strip.background = ggplot2::element_rect(fill = "grey94", color = "grey75", linewidth = 0.35),
      strip.text = ggplot2::element_text(face = "bold", color = "black"),
      legend.key = ggplot2::element_blank(),
      plot.title.position = "plot"
    )
}

format_p <- function(x) {
  ifelse(
    is.na(x), "NA",
    ifelse(x < 0.001, formatC(x, format = "e", digits = 1), sprintf("%.3f", x))
  )
}

cat("Reading Seurat object and frozen pathway definitions...\n")
obj <- readRDS(paths$object)
if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
  obj <- SeuratObject::JoinLayers(obj, assay = "RNA")
}
counts <- Seurat::GetAssayData(obj, assay = "RNA", layer = "counts")
meta <- tibble::rownames_to_column(obj@meta.data, "cell")

cnv <- utils::read.csv(paths$cnv, check.names = FALSE)
cnv <- cnv[, c("cell", "infercnv_mean_abs_log2"), drop = FALSE]
meta <- meta |>
  dplyr::select(-dplyr::any_of("infercnv_mean_abs_log2")) |>
  dplyr::left_join(cnv, by = "cell") |>
  dplyr::mutate(
    patient_pair = factor(patient_pair, levels = c("P1", "P2", "P3")),
    condition_simple = factor(condition_simple, levels = c("diagnosis", "relapse"))
  )

required <- c("cell", "sample_id", "patient_pair", "condition_simple", "broad_annotation",
              "nCount_RNA", "nFeature_RNA", "percent.mt", "infercnv_mean_abs_log2")
if (!all(required %in% colnames(meta))) {
  stop("Missing metadata columns: ", paste(setdiff(required, colnames(meta)), collapse = ", "), call. = FALSE)
}

candidate <- meta |>
  dplyr::filter(
    broad_annotation == "T_lineage_candidate_malignant",
    is.finite(infercnv_mean_abs_log2)
  )

sample_info_all <- candidate |>
  dplyr::distinct(sample_id, patient_pair, condition_simple) |>
  dplyr::arrange(patient_pair, condition_simple)

# Reconstruct the exact frozen KEGG pathway definitions used in the primary run.
kegg_frozen <- utils::read.csv(paths$kegg, check.names = FALSE)
kegg_frozen <- kegg_frozen |>
  dplyr::filter(!duplicated(pathway), nzchar(linked_genes))
kegg_sets <- strsplit(kegg_frozen$linked_genes, ";", fixed = TRUE)
kegg_sets <- lapply(kegg_sets, function(x) sort(unique(trimws(x[nzchar(trimws(x))]))))
names(kegg_sets) <- kegg_frozen$pathway
kegg_annotation <- kegg_frozen |>
  dplyr::select(pathway, pathway_id, pathway_name, gene_set_size_total)

cnv_threshold <- function(top_fraction) {
  unname(stats::quantile(candidate$infercnv_mean_abs_log2, 1 - top_fraction, na.rm = TRUE, type = 7))
}
primary_threshold <- cnv_threshold(0.25)

sample_qc_bounds <- candidate |>
  dplyr::group_by(sample_id) |>
  dplyr::summarise(
    nf_q95 = unname(stats::quantile(nFeature_RNA, 0.95, na.rm = TRUE)),
    nf_q99 = unname(stats::quantile(nFeature_RNA, 0.99, na.rm = TRUE)),
    .groups = "drop"
  )
candidate <- candidate |> dplyr::left_join(sample_qc_bounds, by = "sample_id")

case_specs <- tibble::tribble(
  ~case_id, ~axis, ~case_label, ~selection_mode, ~top_fraction, ~qc_rule, ~rank_metric, ~gene_filter, ~robust, ~norm_method, ~leave_out,
  "primary", "Primary", "Primary: global CNV top 25%", "global", .25, "original", "signed_sqrt_F", "filterByExpr", TRUE, "TMM", "",
  "cnv_top15", "CNV cutoff", "CNV top 15%", "global", .15, "original", "signed_sqrt_F", "filterByExpr", TRUE, "TMM", "",
  "cnv_top20", "CNV cutoff", "CNV top 20%", "global", .20, "original", "signed_sqrt_F", "filterByExpr", TRUE, "TMM", "",
  "cnv_top30", "CNV cutoff", "CNV top 30%", "global", .30, "original", "signed_sqrt_F", "filterByExpr", TRUE, "TMM", "",
  "cnv_top40", "CNV cutoff", "CNV top 40%", "global", .40, "original", "signed_sqrt_F", "filterByExpr", TRUE, "TMM", "",
  "cnv_top50", "CNV cutoff", "CNV top 50%", "global", .50, "original", "signed_sqrt_F", "filterByExpr", TRUE, "TMM", "",
  "cnv_sample_top25", "CNV cutoff", "Top 25% within each sample", "within_sample", .25, "original", "signed_sqrt_F", "filterByExpr", TRUE, "TMM", "",
  "qc_mt4", "QC", "Mitochondrial RNA < 4%", "global", .25, "mt4", "signed_sqrt_F", "filterByExpr", TRUE, "TMM", "",
  "qc_mt3", "QC", "Mitochondrial RNA < 3%", "global", .25, "mt3", "signed_sqrt_F", "filterByExpr", TRUE, "TMM", "",
  "qc_nf500", "QC", "Features >= 500", "global", .25, "nf500", "signed_sqrt_F", "filterByExpr", TRUE, "TMM", "",
  "qc_nf1000", "QC", "Features >= 1,000", "global", .25, "nf1000", "signed_sqrt_F", "filterByExpr", TRUE, "TMM", "",
  "qc_nf1500", "QC", "Features >= 1,500", "global", .25, "nf1500", "signed_sqrt_F", "filterByExpr", TRUE, "TMM", "",
  "qc_upper99", "QC", "Remove top 1% feature counts/sample", "global", .25, "nf_upper99", "signed_sqrt_F", "filterByExpr", TRUE, "TMM", "",
  "qc_upper95", "QC", "Remove top 5% feature counts/sample", "global", .25, "nf_upper95", "signed_sqrt_F", "filterByExpr", TRUE, "TMM", "",
  "qc_combined", "QC", "Combined strict QC", "global", .25, "strict_combo", "signed_sqrt_F", "filterByExpr", TRUE, "TMM", "",
  "filter_cpm05_2", "Gene filter", "CPM >= 0.5 in >= 2 samples", "global", .25, "original", "signed_sqrt_F", "cpm05_n2", TRUE, "TMM", "",
  "filter_cpm1_2", "Gene filter", "CPM >= 1 in >= 2 samples", "global", .25, "original", "signed_sqrt_F", "cpm1_n2", TRUE, "TMM", "",
  "filter_cpm1_3", "Gene filter", "CPM >= 1 in >= 3 samples", "global", .25, "original", "signed_sqrt_F", "cpm1_n3", TRUE, "TMM", "",
  "model_nonrobust", "Model", "Non-robust dispersion", "global", .25, "original", "signed_sqrt_F", "filterByExpr", FALSE, "TMM", "",
  "model_tmmwsp", "Model", "TMMwsp normalization", "global", .25, "original", "signed_sqrt_F", "filterByExpr", TRUE, "TMMwsp", "",
  "rank_logfc", "GSEA ranking", "Rank by signed log2FC", "global", .25, "original", "logFC", "filterByExpr", TRUE, "TMM", "",
  "rank_signed_p", "GSEA ranking", "Rank by signed -log10(P)", "global", .25, "original", "signed_neglog10_P", "filterByExpr", TRUE, "TMM", "",
  "leaveout_p1", "Patient influence", "Leave out patient 1", "global", .25, "original", "signed_sqrt_F", "filterByExpr", TRUE, "TMM", "P1",
  "leaveout_p2", "Patient influence", "Leave out patient 2", "global", .25, "original", "signed_sqrt_F", "filterByExpr", TRUE, "TMM", "P2",
  "leaveout_p3", "Patient influence", "Leave out patient 3", "global", .25, "original", "signed_sqrt_F", "filterByExpr", TRUE, "TMM", "P3"
)

select_case_cells <- function(spec) {
  x <- candidate
  if (identical(spec$selection_mode, "within_sample")) {
    x <- x |>
      dplyr::group_by(sample_id) |>
      dplyr::mutate(cnv_case_threshold = unname(stats::quantile(
        infercnv_mean_abs_log2, 1 - spec$top_fraction, na.rm = TRUE, type = 7
      ))) |>
      dplyr::ungroup() |>
      dplyr::filter(infercnv_mean_abs_log2 >= cnv_case_threshold)
  } else {
    threshold <- cnv_threshold(spec$top_fraction)
    x <- x |> dplyr::filter(infercnv_mean_abs_log2 >= threshold)
  }

  x <- switch(
    spec$qc_rule,
    original = x,
    mt4 = dplyr::filter(x, percent.mt < 4),
    mt3 = dplyr::filter(x, percent.mt < 3),
    nf500 = dplyr::filter(x, nFeature_RNA >= 500),
    nf1000 = dplyr::filter(x, nFeature_RNA >= 1000),
    nf1500 = dplyr::filter(x, nFeature_RNA >= 1500),
    nf_upper99 = dplyr::filter(x, nFeature_RNA <= nf_q99),
    nf_upper95 = dplyr::filter(x, nFeature_RNA <= nf_q95),
    strict_combo = dplyr::filter(x, percent.mt < 4, nFeature_RNA >= 500, nFeature_RNA <= nf_q99),
    stop("Unknown QC rule: ", spec$qc_rule)
  )
  if (nzchar(spec$leave_out)) x <- dplyr::filter(x, as.character(patient_pair) != spec$leave_out)
  x
}

make_pb <- function(meta_case, sample_info) {
  cells <- meta_case$cell
  sample_index <- match(meta_case$sample_id, sample_info$sample_id)
  map <- Matrix::sparseMatrix(
    i = seq_along(cells), j = sample_index, x = 1,
    dims = c(length(cells), nrow(sample_info)),
    dimnames = list(cells, sample_info$sample_id)
  )
  as(counts[, cells, drop = FALSE] %*% map, "dgCMatrix")
}

filter_genes <- function(y, design, filter_name) {
  if (filter_name == "filterByExpr") return(edgeR::filterByExpr(y, design = design))
  cpm_mat <- edgeR::cpm(y)
  if (filter_name == "cpm05_n2") return(Matrix::rowSums(cpm_mat >= 0.5) >= 2)
  if (filter_name == "cpm1_n2") return(Matrix::rowSums(cpm_mat >= 1) >= 2)
  if (filter_name == "cpm1_n3") return(Matrix::rowSums(cpm_mat >= 1) >= 3)
  stop("Unknown gene filter: ", filter_name)
}

rank_from_deg <- function(deg, metric) {
  out <- switch(
    metric,
    signed_sqrt_F = sign(deg$logFC) * sqrt(pmax(deg$F, 0)),
    logFC = deg$logFC,
    signed_neglog10_P = sign(deg$logFC) * -log10(pmax(deg$PValue, 1e-300)),
    stop("Unknown ranking metric: ", metric)
  )
  names(out) <- deg$gene
  out <- out[is.finite(out) & !duplicated(names(out))]
  # Deterministic epsilon breaks exact ties without changing meaningful order.
  out <- out + rank(names(out), ties.method = "first") * 1e-12
  sort(out, decreasing = TRUE)
}

run_case <- function(spec) {
  meta_case <- select_case_cells(spec)
  sample_info <- sample_info_all
  if (nzchar(spec$leave_out)) {
    sample_info <- sample_info |> dplyr::filter(as.character(patient_pair) != spec$leave_out)
    sample_info$patient_pair <- droplevels(sample_info$patient_pair)
  }
  sample_info$condition_simple <- droplevels(sample_info$condition_simple)

  cell_counts <- sample_info |>
    dplyr::left_join(meta_case |> dplyr::count(sample_id, name = "n_cells"), by = "sample_id") |>
    dplyr::mutate(
      n_cells = ifelse(is.na(n_cells), 0L, n_cells),
      case_id = spec$case_id,
      axis = spec$axis,
      case_label = spec$case_label,
      .before = 1
    )
  valid <- nrow(sample_info) >= 4 && all(cell_counts$n_cells >= 20)
  registry <- spec |>
    dplyr::mutate(
      n_patients = dplyr::n_distinct(sample_info$patient_pair),
      n_samples = nrow(sample_info),
      n_cells_total = nrow(meta_case),
      min_cells_per_sample = min(cell_counts$n_cells),
      median_cells_per_sample = stats::median(cell_counts$n_cells),
      valid = valid,
      invalid_reason = ifelse(valid, "", "requires at least 20 cells in every included sample")
    )
  if (!valid) return(list(registry = registry, cell_counts = cell_counts, result = NULL))

  pb <- make_pb(meta_case, sample_info)
  y <- edgeR::DGEList(counts = pb, samples = as.data.frame(sample_info))
  design <- stats::model.matrix(~ patient_pair + condition_simple, data = sample_info)
  keep <- filter_genes(y, design, spec$gene_filter)
  if (sum(keep) < 100) {
    registry$valid <- FALSE
    registry$invalid_reason <- "fewer than 100 genes after filtering"
    return(list(registry = registry, cell_counts = cell_counts, result = NULL))
  }
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y, method = spec$norm_method)
  y <- edgeR::estimateDisp(y, design, robust = spec$robust)
  fit <- edgeR::glmQLFit(y, design, robust = spec$robust)
  coef_name <- grep("condition_simple.*relapse", colnames(design), value = TRUE)
  if (length(coef_name) != 1) stop("Could not identify relapse coefficient for ", spec$case_id)
  qlf <- edgeR::glmQLFTest(fit, coef = coef_name)
  deg <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
  deg$gene <- rownames(deg)
  deg$FDR <- stats::p.adjust(deg$PValue, method = "BH")
  rank_stats <- rank_from_deg(deg, spec$rank_metric)

  tested_sets <- lapply(kegg_sets, intersect, y = names(rank_stats))
  tested_sets <- tested_sets[vapply(tested_sets, length, integer(1)) >= 5]
  set.seed(20260909 + match(spec$case_id, case_specs$case_id))
  fg <- suppressWarnings(fgsea::fgseaMultilevel(
    pathways = tested_sets, stats = rank_stats,
    minSize = 5, maxSize = 500, eps = 0, nproc = 1
  ))
  fg <- as.data.frame(fg)
  fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  fg <- fg |>
    dplyr::left_join(kegg_annotation, by = "pathway") |>
    dplyr::mutate(
      case_id = spec$case_id,
      axis = spec$axis,
      case_label = spec$case_label,
      primary_analysis = spec$case_id == "primary",
      n_patients = dplyr::n_distinct(sample_info$patient_pair),
      n_cells_total = nrow(meta_case),
      min_cells_per_sample = min(cell_counts$n_cells),
      genes_ranked_n = length(rank_stats),
      direction = dplyr::case_when(NES > 0 ~ "relapse_up", NES < 0 ~ "relapse_down", TRUE ~ "flat"),
      global_FDR_tier = dplyr::case_when(
        padj < 0.05 ~ "FDR < 0.05",
        padj < 0.10 ~ "FDR 0.05-0.10",
        padj < 0.20 ~ "FDR 0.10-0.20",
        TRUE ~ "FDR >= 0.20"
      ),
      .before = 1
    )

  ids2 <- c("hsa03050", "hsa04137")
  ids4 <- c(ids2, "hsa04136", "hsa04140")
  fg$target_family2_BH <- NA_real_
  fg$target_family2_Holm <- NA_real_
  fg$target_family4_BH <- NA_real_
  fg$target_family4_Holm <- NA_real_
  i2 <- which(fg$pathway_id %in% ids2)
  i4 <- which(fg$pathway_id %in% ids4)
  fg$target_family2_BH[i2] <- stats::p.adjust(fg$pval[i2], method = "BH")
  fg$target_family2_Holm[i2] <- stats::p.adjust(fg$pval[i2], method = "holm")
  fg$target_family4_BH[i4] <- stats::p.adjust(fg$pval[i4], method = "BH")
  fg$target_family4_Holm[i4] <- stats::p.adjust(fg$pval[i4], method = "holm")

  list(
    registry = registry,
    cell_counts = cell_counts,
    result = fg,
    deg = if (spec$case_id == "primary") deg else NULL,
    logcpm = if (spec$case_id == "primary") edgeR::cpm(y, log = TRUE, prior.count = 2) else NULL,
    sample_info = if (spec$case_id == "primary") sample_info else NULL
  )
}

reuse_existing <- tolower(get_arg("--reuse-results", "false")) %in% c("true", "1", "yes")
registry_file <- file.path(paths$tables, "robustness_parameter_registry.csv")
cell_counts_file <- file.path(paths$tables, "robustness_cell_counts_by_sample.csv")
all_gsea_file <- file.path(paths$tables, "robustness_KEGG_GSEA_all_cases.csv")

if (reuse_existing && all(file.exists(c(registry_file, cell_counts_file, all_gsea_file)))) {
  cat("Reusing completed sensitivity grid; recomputing primary pseudobulk for paired score figure...\n")
  registry <- utils::read.csv(registry_file, check.names = FALSE)
  cell_counts <- utils::read.csv(cell_counts_file, check.names = FALSE)
  all_gsea <- utils::read.csv(all_gsea_file, check.names = FALSE)
  registry$axis <- case_specs$axis[match(registry$case_id, case_specs$case_id)]
  registry$case_label <- case_specs$case_label[match(registry$case_id, case_specs$case_id)]
  cell_counts$axis <- case_specs$axis[match(cell_counts$case_id, case_specs$case_id)]
  cell_counts$case_label <- case_specs$case_label[match(cell_counts$case_id, case_specs$case_id)]
  all_gsea$axis <- case_specs$axis[match(all_gsea$case_id, case_specs$case_id)]
  all_gsea$case_label <- case_specs$case_label[match(all_gsea$case_id, case_specs$case_id)]
  primary_run <- run_case(case_specs[case_specs$case_id == "primary", , drop = FALSE])
} else {
  cat("Running", nrow(case_specs), "predefined sensitivity cases...\n")
  case_results <- vector("list", nrow(case_specs))
  for (i in seq_len(nrow(case_specs))) {
    cat(sprintf("[%02d/%02d] %s\n", i, nrow(case_specs), case_specs$case_id[[i]]))
    case_results[[i]] <- run_case(case_specs[i, , drop = FALSE])
  }
  registry <- dplyr::bind_rows(lapply(case_results, `[[`, "registry"))
  cell_counts <- dplyr::bind_rows(lapply(case_results, `[[`, "cell_counts"))
  all_gsea <- dplyr::bind_rows(lapply(case_results, `[[`, "result"))
  primary_run <- case_results[[match("primary", case_specs$case_id)]]
}

# Preserve the originally frozen primary GSEA result. A fresh fgseaMultilevel
# call is stochastic near the threshold, so its single-seed recomputation is
# stored as a Monte Carlo audit rather than allowed to replace the main result.
primary_recomputed <- primary_run$result |>
  dplyr::select(pathway_id, recomputed_NES = NES, recomputed_nominal_P = pval, recomputed_global_FDR = padj)
mc_compare <- kegg_frozen |>
  dplyr::filter(pathway_id %in% c("hsa03050", "hsa04137", "hsa04136", "hsa04140")) |>
  dplyr::select(pathway_id, pathway_name, frozen_NES = NES, frozen_nominal_P = pval, frozen_global_FDR = padj) |>
  dplyr::left_join(primary_recomputed, by = "pathway_id") |>
  dplyr::mutate(
    delta_NES = recomputed_NES - frozen_NES,
    nominal_P_ratio_recomputed_to_frozen = recomputed_nominal_P / frozen_nominal_P,
    note = "Single-seed fgseaMultilevel recomputation; the frozen primary result remains authoritative."
  )
write_csv_file(mc_compare, file.path(paths$tables, "fgsea_primary_recompute_monte_carlo_check.csv"))

primary_index_all <- which(all_gsea$case_id == "primary")
frozen_match <- match(all_gsea$pathway[primary_index_all], kegg_frozen$pathway)
keep_frozen <- !is.na(frozen_match)
primary_rows <- primary_index_all[keep_frozen]
frozen_match <- frozen_match[keep_frozen]
for (nm in intersect(c("pval", "padj", "ES", "NES", "size", "leadingEdge", "log2err"), colnames(all_gsea))) {
  all_gsea[primary_rows, nm] <- kegg_frozen[frozen_match, nm]
}
all_gsea$direction[primary_rows] <- ifelse(all_gsea$NES[primary_rows] > 0, "relapse_up", "relapse_down")
all_gsea$global_FDR_tier[primary_rows] <- dplyr::case_when(
  all_gsea$padj[primary_rows] < 0.05 ~ "FDR < 0.05",
  all_gsea$padj[primary_rows] < 0.10 ~ "FDR 0.05-0.10",
  all_gsea$padj[primary_rows] < 0.20 ~ "FDR 0.10-0.20",
  TRUE ~ "FDR >= 0.20"
)
for (nm in c("target_family2_BH", "target_family2_Holm", "target_family4_BH", "target_family4_Holm")) {
  all_gsea[primary_rows, nm] <- NA_real_
}
i2 <- which(all_gsea$case_id == "primary" & all_gsea$pathway_id %in% c("hsa03050", "hsa04137"))
i4 <- which(all_gsea$case_id == "primary" & all_gsea$pathway_id %in% c("hsa03050", "hsa04137", "hsa04136", "hsa04140"))
all_gsea$target_family2_BH[i2] <- stats::p.adjust(all_gsea$pval[i2], method = "BH")
all_gsea$target_family2_Holm[i2] <- stats::p.adjust(all_gsea$pval[i2], method = "holm")
all_gsea$target_family4_BH[i4] <- stats::p.adjust(all_gsea$pval[i4], method = "BH")
all_gsea$target_family4_Holm[i4] <- stats::p.adjust(all_gsea$pval[i4], method = "holm")

focus_ids <- c("hsa05012", "hsa03010", "hsa00190", "hsa04630", "hsa03050", "hsa04137", "hsa04136", "hsa04140")
focus_gsea <- all_gsea |>
  dplyr::filter(pathway_id %in% focus_ids) |>
  dplyr::arrange(match(case_id, case_specs$case_id), match(pathway_id, focus_ids))

write_csv_file(registry, registry_file)
write_csv_file(cell_counts, cell_counts_file)
write_csv_file(all_gsea, all_gsea_file)
write_csv_file(focus_gsea, file.path(paths$tables, "robustness_KEGG_GSEA_focus_pathways.csv"))

# Selection diagnostics explain why P values change when definitions or QC move.
global_thresholds <- tibble::tibble(top_fraction = c(.15, .20, .25, .30, .40, .50)) |>
  dplyr::mutate(
    threshold_scope = "global candidate pool",
    sample_id = "all",
    cnv_threshold = vapply(top_fraction, cnv_threshold, numeric(1)),
    candidate_cells_n = nrow(candidate)
  )
sample_thresholds <- candidate |>
  dplyr::group_by(sample_id, patient_pair, condition_simple) |>
  dplyr::summarise(
    threshold_scope = "within sample",
    top_fraction = .25,
    cnv_threshold = unname(stats::quantile(infercnv_mean_abs_log2, .75, na.rm = TRUE)),
    candidate_cells_n = dplyr::n(),
    .groups = "drop"
  )
write_csv_file(
  dplyr::bind_rows(global_thresholds, sample_thresholds),
  file.path(paths$tables, "robustness_CNV_selection_thresholds.csv")
)

selection_by_sample <- dplyr::bind_rows(lapply(seq_len(nrow(case_specs)), function(i) {
  spec <- case_specs[i, , drop = FALSE]
  x <- select_case_cells(spec)
  sample_info <- sample_info_all
  if (nzchar(spec$leave_out)) sample_info <- sample_info |> dplyr::filter(as.character(patient_pair) != spec$leave_out)
  sample_info |>
    dplyr::left_join(
      x |>
        dplyr::group_by(sample_id) |>
        dplyr::summarise(
          n_cells = dplyr::n(),
          median_nFeature_RNA = stats::median(nFeature_RNA),
          median_percent_mt = stats::median(percent.mt),
          median_CNV_burden = stats::median(infercnv_mean_abs_log2),
          .groups = "drop"
        ),
      by = "sample_id"
    ) |>
    dplyr::mutate(case_id = spec$case_id, axis = spec$axis, case_label = spec$case_label, .before = 1)
}))
write_csv_file(selection_by_sample, file.path(paths$tables, "robustness_selection_QC_by_sample.csv"))

state_composition <- dplyr::bind_rows(lapply(seq_len(nrow(case_specs)), function(i) {
  spec <- case_specs[i, , drop = FALSE]
  x <- select_case_cells(spec)
  x |>
    dplyr::count(major_annotation, name = "n_cells") |>
    dplyr::mutate(
      proportion = n_cells / sum(n_cells),
      case_id = spec$case_id,
      axis = spec$axis,
      case_label = spec$case_label,
      .before = 1
    )
}))
write_csv_file(state_composition, file.path(paths$tables, "robustness_cell_state_composition.csv"))

# Multiplicity table: global KEGG FDR is the primary exploratory-screen result.
# Target-family corrections are illustrative and only confirmatory if the family
# was fixed before reviewing these data.
primary_mult <- kegg_frozen |>
  dplyr::filter(pathway_id %in% c("hsa03050", "hsa04137", "hsa04136", "hsa04140")) |>
  dplyr::transmute(
    pathway_id, pathway_name, NES, nominal_P = pval,
    global_KEGG_FDR = padj
  )
primary_mult$BH_two_target_family <- NA_real_
primary_mult$Holm_two_target_family <- NA_real_
primary_mult$BH_four_pathway_family <- NA_real_
primary_mult$Holm_four_pathway_family <- NA_real_
i2 <- which(primary_mult$pathway_id %in% c("hsa03050", "hsa04137"))
i4 <- which(primary_mult$pathway_id %in% c("hsa03050", "hsa04137", "hsa04136", "hsa04140"))
primary_mult$BH_two_target_family[i2] <- stats::p.adjust(primary_mult$nominal_P[i2], method = "BH")
primary_mult$Holm_two_target_family[i2] <- stats::p.adjust(primary_mult$nominal_P[i2], method = "holm")
primary_mult$BH_four_pathway_family[i4] <- stats::p.adjust(primary_mult$nominal_P[i4], method = "BH")
primary_mult$Holm_four_pathway_family[i4] <- stats::p.adjust(primary_mult$nominal_P[i4], method = "holm")
primary_mult <- primary_mult |>
  dplyr::mutate(
    interpretation = dplyr::case_when(
      global_KEGG_FDR < 0.05 ~ "global KEGG FDR < 0.05",
      global_KEGG_FDR < 0.10 ~ "exploratory global KEGG FDR < 0.10",
      TRUE ~ "not global KEGG FDR < 0.10"
    ),
    preregistration_caveat = "Target-family correction is confirmatory only when hypotheses and analysis plan were fixed before data review."
  )
write_csv_file(primary_mult, file.path(paths$tables, "primary_GSEA_multiple_testing_framework.csv"))

# Primary leading-edge patient-direction audit.
primary_path <- kegg_frozen |>
  dplyr::filter(pathway_id %in% c("hsa03050", "hsa04137", "hsa04136", "hsa04140"))
pair_score_rows <- lapply(seq_len(nrow(primary_path)), function(i) {
  genes <- strsplit(primary_path$leadingEdge[[i]], ";", fixed = TRUE)[[1]]
  genes <- intersect(genes, rownames(primary_run$logcpm))
  if (!length(genes)) return(NULL)
  score <- Matrix::colMeans(primary_run$logcpm[genes, , drop = FALSE])
  tibble::tibble(
    pathway_id = primary_path$pathway_id[[i]],
    pathway_name = primary_path$pathway_name[[i]],
    GSEA_NES = primary_path$NES[[i]],
    GSEA_nominal_P = primary_path$pval[[i]],
    GSEA_global_FDR = primary_path$padj[[i]],
    leading_edge_genes_n = length(genes),
    sample_id = names(score),
    leading_edge_mean_logCPM = as.numeric(score)
  ) |>
    dplyr::left_join(primary_run$sample_info, by = "sample_id")
})
pair_scores <- dplyr::bind_rows(pair_score_rows) |>
  dplyr::group_by(pathway_id, pathway_name) |>
  dplyr::mutate(
    paired_delta = leading_edge_mean_logCPM[condition_simple == "relapse"][match(patient_pair, patient_pair[condition_simple == "relapse"])] -
      leading_edge_mean_logCPM[condition_simple == "diagnosis"][match(patient_pair, patient_pair[condition_simple == "diagnosis"])]
  ) |>
  dplyr::ungroup()
pair_summary <- pair_scores |>
  dplyr::distinct(pathway_id, pathway_name, patient_pair, paired_delta, GSEA_NES, GSEA_nominal_P, GSEA_global_FDR) |>
  dplyr::group_by(pathway_id, pathway_name, GSEA_NES, GSEA_nominal_P, GSEA_global_FDR) |>
  dplyr::summarise(
    n_pairs = dplyr::n(),
    pairs_relapse_higher = sum(paired_delta > 0),
    median_paired_delta = stats::median(paired_delta),
    paired_t_P = suppressWarnings(stats::t.test(paired_delta, mu = 0)$p.value),
    exact_paired_Wilcoxon_P = suppressWarnings(stats::wilcox.test(paired_delta, mu = 0, paired = FALSE, exact = TRUE)$p.value),
    .groups = "drop"
  )
write_csv_file(pair_scores, file.path(paths$tables, "primary_leading_edge_scores_by_sample.csv"))
write_csv_file(pair_summary, file.path(paths$tables, "primary_leading_edge_pair_direction_summary.csv"))

# Robustness summary across every predefined, valid case. This is descriptive;
# it does not turn the smallest observed P value into the primary result.
robust_summary <- focus_gsea |>
  dplyr::filter(pathway_id %in% c("hsa03050", "hsa04137", "hsa04136", "hsa04140")) |>
  dplyr::group_by(pathway_id, pathway_name) |>
  dplyr::summarise(
    valid_cases_n = dplyr::n(),
    positive_NES_cases_n = sum(NES > 0),
    negative_NES_cases_n = sum(NES < 0),
    nominal_P_lt_0.05_cases_n = sum(pval < 0.05),
    global_FDR_lt_0.10_cases_n = sum(padj < 0.10),
    global_FDR_lt_0.05_cases_n = sum(padj < 0.05),
    NES_min = min(NES), NES_median = stats::median(NES), NES_max = max(NES),
    nominal_P_min = min(pval), nominal_P_max = max(pval),
    global_FDR_min = min(padj), global_FDR_max = max(padj),
    .groups = "drop"
  )
write_csv_file(robust_summary, file.path(paths$tables, "robustness_focus_pathway_summary.csv"))

# Same-direction cross-dataset ranking uses the worse of the two FDR values so
# a pathway cannot rank highly because it is strong in only one dataset.
if (file.exists(paths$cross)) {
  cross <- utils::read.csv(paths$cross, check.names = FALSE) |>
    dplyr::mutate(
      worst_dataset_FDR = pmax(sc_FDR, bulk_FDR),
      geometric_mean_FDR = sqrt(sc_FDR * bulk_FDR),
      mean_NES = (sc_NES + bulk_NES) / 2
    ) |>
    dplyr::arrange(worst_dataset_FDR, geometric_mean_FDR) |>
    dplyr::mutate(robust_shared_rank = dplyr::row_number())
  write_csv_file(cross, file.path(paths$tables, "cross_dataset_same_direction_pathway_robustness_rank.csv"))
}

# Figure 1: all sensitivity cases for the two central hypotheses.
plot_dat <- focus_gsea |>
  dplyr::filter(pathway_id %in% c("hsa03050", "hsa04137")) |>
  dplyr::mutate(
    pathway_name = factor(pathway_name, levels = c("Mitophagy - animal", "Proteasome")),
    axis = factor(axis, levels = unique(case_specs$axis)),
    case_label = factor(case_label, levels = rev(case_specs$case_label)),
    evidence_size = pmin(-log10(pmax(pval, 1e-12)), 5),
    global_FDR_tier = dplyr::case_when(
      padj < 0.05 ~ "FDR < 0.05",
      padj < 0.10 ~ "FDR 0.05-0.10",
      padj < 0.20 ~ "FDR 0.10-0.20",
      TRUE ~ "FDR >= 0.20"
    ),
    global_FDR_tier = factor(global_FDR_tier, levels = c("FDR < 0.05", "FDR 0.05-0.10", "FDR 0.10-0.20", "FDR >= 0.20"))
  )

p_robust <- ggplot2::ggplot(plot_dat, ggplot2::aes(
  x = NES, y = case_label, color = pathway_name,
  shape = global_FDR_tier, size = evidence_size
)) +
  ggplot2::geom_vline(xintercept = 0, color = "grey75", linewidth = 0.4) +
  ggplot2::geom_point(alpha = 0.95, stroke = 0.9) +
  ggplot2::facet_grid(axis ~ ., scales = "free_y", space = "free_y", switch = "y") +
  ggplot2::scale_color_manual(values = c("Mitophagy - animal" = "#0072B2", "Proteasome" = "#D55E00")) +
  ggplot2::scale_shape_manual(values = c("FDR < 0.05" = 16, "FDR 0.05-0.10" = 17, "FDR 0.10-0.20" = 15, "FDR >= 0.20" = 1)) +
  ggplot2::scale_size_continuous(name = expression(-log[10]("nominal P")), range = c(2.2, 5.2), breaks = c(1, 2, 3, 4)) +
  ggplot2::labs(
    title = "Proteasome and mitophagy across predefined sensitivity analyses",
    subtitle = "Paired pseudobulk; every valid setting is shown. Shape encodes global KEGG FDR tier.",
    x = "Normalized enrichment score (relapse vs diagnosis)", y = NULL,
    color = "KEGG pathway", shape = "Global KEGG FDR"
  ) +
  theme_pub(9) +
  ggplot2::theme(
    strip.placement = "outside",
    strip.text.y.left = ggplot2::element_text(angle = 0),
    axis.text.y = ggplot2::element_text(size = 8.4),
    legend.position = "bottom",
    legend.box = "vertical",
    panel.spacing.y = grid::unit(0.55, "lines"),
    plot.margin = ggplot2::margin(8, 12, 8, 8)
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(order = 1, nrow = 1),
    shape = ggplot2::guide_legend(order = 2, nrow = 1),
    size = ggplot2::guide_legend(order = 3, nrow = 1)
  )
ggplot2::ggsave(
  file.path(paths$figures, "FigSens04_KEGG_GSEA_parameter_robustness.png"),
  p_robust, width = 9.2, height = 13.2, dpi = 360, bg = "white"
)

# Figure 2: the same primary P values under different multiplicity families.
mult_long <- primary_mult |>
  dplyr::select(pathway_name, nominal_P, BH_two_target_family, BH_four_pathway_family, global_KEGG_FDR) |>
  tidyr::pivot_longer(-pathway_name, names_to = "correction", values_to = "P_or_FDR") |>
  dplyr::mutate(
    correction = factor(
      correction,
      levels = c("nominal_P", "BH_two_target_family", "BH_four_pathway_family", "global_KEGG_FDR"),
      labels = c("Nominal P", "BH: 2 targets", "BH: 4 pathways", "Global KEGG FDR")
    ),
    pathway_name = factor(pathway_name, levels = rev(c("Proteasome", "Mitophagy - animal", "Autophagy - other", "Autophagy - animal"))),
    evidence = pmin(-log10(pmax(P_or_FDR, 1e-12)), 5),
    label = format_p(P_or_FDR)
  )
p_mult <- ggplot2::ggplot(mult_long, ggplot2::aes(correction, pathway_name, fill = evidence)) +
  ggplot2::geom_tile(color = "white", linewidth = 1) +
  ggplot2::geom_text(ggplot2::aes(label = label), size = 3.4, color = "black") +
  ggplot2::scale_fill_gradient(low = "#F2F2F2", high = "#B2182B", na.value = "grey88", name = expression(-log[10](P))) +
  ggplot2::labs(
    title = "Multiplicity changes the strength of the claim, not the effect direction",
    subtitle = "Target-family columns are confirmatory only if fixed before examining these data.",
    x = NULL, y = NULL
  ) +
  theme_pub(10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 18, hjust = 1),
    legend.position = "right",
    axis.line = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank()
  )
ggplot2::ggsave(
  file.path(paths$figures, "FigSens05_GSEA_multiple_testing_framework.png"),
  p_mult, width = 8.4, height = 4.6, dpi = 360, bg = "white"
)

# Figure 3: paired patient-level direction of primary leading-edge scores.
pair_annot <- pair_summary |>
  dplyr::mutate(
    label = paste0("Relapse higher: ", pairs_relapse_higher, "/", n_pairs,
                   "\nGSEA FDR = ", format_p(GSEA_global_FDR)),
    condition_simple = "relapse"
  )
p_pair <- ggplot2::ggplot(pair_scores, ggplot2::aes(
  x = condition_simple, y = leading_edge_mean_logCPM, group = patient_pair
)) +
  ggplot2::geom_line(color = "grey60", linewidth = 0.55) +
  ggplot2::geom_point(
    ggplot2::aes(fill = condition_simple, shape = condition_simple),
    color = "black", size = 3.0, stroke = 0.55
  ) +
  ggplot2::geom_text(
    data = pair_annot,
    ggplot2::aes(x = condition_simple, y = Inf, label = label),
    inherit.aes = FALSE, hjust = 1.05, vjust = 1.25, size = 3.0, lineheight = 0.95
  ) +
  ggplot2::facet_wrap(~ pathway_name, scales = "free_y", ncol = 2) +
  ggplot2::scale_x_discrete(limits = c("diagnosis", "relapse"), labels = c("Diagnosis", "Relapse")) +
  ggplot2::scale_fill_manual(values = c("diagnosis" = "#56B4E9", "relapse" = "#D55E00")) +
  ggplot2::scale_shape_manual(values = c("diagnosis" = 21, "relapse" = 24)) +
  ggplot2::labs(
    title = "Patient-paired direction of primary GSEA leading-edge programs",
    subtitle = "Each line is one patient; annotation reports pathway-level GSEA FDR, not a cell-level test.",
    x = NULL, y = "Mean pseudobulk logCPM across leading-edge genes"
  ) +
  theme_pub(10) +
  ggplot2::theme(legend.position = "none", panel.spacing = grid::unit(1.0, "lines"))
ggplot2::ggsave(
  file.path(paths$figures, "FigSens06_primary_leading_edge_paired_direction.png"),
  p_pair, width = 8.5, height = 7.0, dpi = 360, bg = "white"
)

capture.output(
  sessionInfo(),
  file = file.path(paths$environment, "sessionInfo_17_robustness_sensitivity_audit.txt")
)
cat("Robustness audit complete.\n")
