options(stringsAsFactors = FALSE)
options(timeout = 300)

# Full-KEGG audit of a small, pre-specified set of publication-defensible
# population/QC definitions. All candidates retain the paired patient-level
# pseudobulk model, filterByExpr, robust edgeR QL fitting, TMM normalization,
# and sign(logFC) * sqrt(F) ranking. FDR is BH-adjusted across all tested KEGG
# pathways within each candidate analysis.

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
  tables = file.path(project, "02_Tables", "03_Sensitivity_Checks")
)
stopifnot(file.exists(paths$object), file.exists(paths$cnv), file.exists(paths$kegg))
dir.create(paths$tables, recursive = TRUE, showWarnings = FALSE)

user_lib <- Sys.getenv("R_LIBS_USER", unset = "C:/Users/User/Documents/R/win-library/4.4")
if (dir.exists(user_lib)) .libPaths(unique(c(normalizePath(user_lib, winslash = "/"), .libPaths())))

pkgs <- c("Seurat", "SeuratObject", "Matrix", "edgeR", "fgsea", "dplyr", "tibble")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))

write_csv_file <- function(x, path) {
  utils::write.csv(as.data.frame(x), path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
}

cat("Reading object and pathway definitions...\n")
obj <- readRDS(paths$object)
if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
  obj <- SeuratObject::JoinLayers(obj, assay = "RNA")
}
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
  dplyr::filter(
    broad_annotation == "T_lineage_candidate_malignant",
    is.finite(infercnv_mean_abs_log2)
  )
sample_info <- candidate |>
  dplyr::distinct(sample_id, patient_pair, condition_simple) |>
  dplyr::arrange(patient_pair, condition_simple)

kegg_frozen <- utils::read.csv(paths$kegg, check.names = FALSE) |>
  dplyr::filter(!duplicated(pathway), nzchar(linked_genes))
kegg_sets <- strsplit(kegg_frozen$linked_genes, ";", fixed = TRUE)
kegg_sets <- lapply(kegg_sets, function(x) sort(unique(trimws(x[nzchar(trimws(x))]))))
names(kegg_sets) <- kegg_frozen$pathway
kegg_annotation <- kegg_frozen |>
  dplyr::select(pathway, pathway_id, pathway_name, gene_set_size_total)

# "none" means the source-paper post-QC object is used unchanged: HTO singlets,
# percent.mt < 5%, and the paper-reported sample-specific lower feature cutoffs.
specs <- tibble::tribble(
  ~candidate_id, ~candidate_label, ~selection_scope, ~top_fraction, ~extra_qc, ~paper_qc_preserved, ~definition_changes,
  "A01", "Global CNV top 25%; paper QC",              "global",        .25, "none",      TRUE,  0L,
  "A02", "Within-sample CNV top 25%; paper QC",       "within_sample", .25, "none",      TRUE,  1L,
  "A03", "Global CNV top 15%; paper QC",              "global",        .15, "none",      TRUE,  1L,
  "A04", "Within-sample CNV top 20%; paper QC",       "within_sample", .20, "none",      TRUE,  2L,
  "A05", "Within-sample CNV top 15%; paper QC",       "within_sample", .15, "none",      TRUE,  2L,
  "A11", "Within-sample CNV top 10%; paper QC",       "within_sample", .10, "none",      TRUE,  2L,
  "A12", "Within-sample CNV top 30%; paper QC",       "within_sample", .30, "none",      TRUE,  2L,
  "A13", "Within-sample CNV top 40%; paper QC",       "within_sample", .40, "none",      TRUE,  2L,
  "A14", "Within-sample CNV top 50%; paper QC",       "within_sample", .50, "none",      TRUE,  2L,
  "A06", "Within-sample CNV top 25%; nFeature >=500", "within_sample", .25, "nf500",     FALSE, 2L,
  "A07", "Within-sample CNV top 15%; nFeature >=500", "within_sample", .15, "nf500",     FALSE, 3L,
  "A08", "Within-sample CNV top 25%; mt <4%",         "within_sample", .25, "mt4",       FALSE, 2L,
  "A09", "Within-sample CNV top 15%; mt <4%",         "within_sample", .15, "mt4",       FALSE, 3L,
  "A10", "Within-sample CNV top 15%; mt <4%; nFeature >=500", "within_sample", .15, "mt4_nf500", FALSE, 4L
)

select_cells <- function(spec) {
  x <- candidate
  if (identical(spec$selection_scope, "within_sample")) {
    x <- x |>
      dplyr::group_by(sample_id) |>
      dplyr::mutate(cnv_threshold = unname(stats::quantile(
        infercnv_mean_abs_log2, 1 - spec$top_fraction, na.rm = TRUE, type = 7
      ))) |>
      dplyr::ungroup() |>
      dplyr::filter(infercnv_mean_abs_log2 >= cnv_threshold)
  } else {
    threshold <- unname(stats::quantile(
      candidate$infercnv_mean_abs_log2, 1 - spec$top_fraction, na.rm = TRUE, type = 7
    ))
    x <- x |> dplyr::filter(infercnv_mean_abs_log2 >= threshold)
  }

  switch(
    spec$extra_qc,
    none = x,
    nf500 = dplyr::filter(x, nFeature_RNA >= 500),
    mt4 = dplyr::filter(x, percent.mt < 4),
    mt4_nf500 = dplyr::filter(x, percent.mt < 4, nFeature_RNA >= 500),
    stop("Unknown extra_qc: ", spec$extra_qc)
  )
}

make_pb <- function(meta_case) {
  cells <- meta_case$cell
  sample_index <- match(meta_case$sample_id, sample_info$sample_id)
  map <- Matrix::sparseMatrix(
    i = seq_along(cells), j = sample_index, x = 1,
    dims = c(length(cells), nrow(sample_info)),
    dimnames = list(cells, sample_info$sample_id)
  )
  as(counts[, cells, drop = FALSE] %*% map, "dgCMatrix")
}

make_rank <- function(deg) {
  out <- sign(deg$logFC) * sqrt(pmax(deg$F, 0))
  names(out) <- deg$gene
  out <- out[is.finite(out) & !duplicated(names(out))]
  out <- out + rank(names(out), ties.method = "first") * 1e-12
  sort(out, decreasing = TRUE)
}

run_candidate <- function(spec, index) {
  cells_meta <- select_cells(spec)
  cell_counts <- sample_info |>
    dplyr::left_join(cells_meta |> dplyr::count(sample_id, name = "n_cells"), by = "sample_id") |>
    dplyr::mutate(n_cells = ifelse(is.na(n_cells), 0L, n_cells))
  if (any(cell_counts$n_cells < 20)) {
    stop(spec$candidate_id, " has fewer than 20 cells in at least one sample")
  }

  pb <- make_pb(cells_meta)
  y <- edgeR::DGEList(counts = pb, samples = as.data.frame(sample_info))
  design <- stats::model.matrix(~ patient_pair + condition_simple, data = sample_info)
  keep <- edgeR::filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y, method = "TMM")
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  coef_name <- grep("condition_simple.*relapse", colnames(design), value = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = coef_name)
  deg <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
  deg$gene <- rownames(deg)
  rank_stats <- make_rank(deg)

  tested_sets <- lapply(kegg_sets, intersect, y = names(rank_stats))
  tested_sets <- tested_sets[vapply(tested_sets, length, integer(1)) >= 5]
  set.seed(20260920L + index)
  fg <- suppressWarnings(fgsea::fgseaMultilevel(
    pathways = tested_sets,
    stats = rank_stats,
    minSize = 5,
    maxSize = 500,
    eps = 0,
    nproc = 1
  )) |>
    as.data.frame()
  fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  result <- fg |>
    dplyr::left_join(kegg_annotation, by = "pathway") |>
    dplyr::mutate(
      candidate_id = spec$candidate_id,
      candidate_label = spec$candidate_label,
      selection_scope = spec$selection_scope,
      top_fraction = spec$top_fraction,
      extra_qc = spec$extra_qc,
      paper_qc_preserved = spec$paper_qc_preserved,
      definition_changes = spec$definition_changes,
      n_cells_total = nrow(cells_meta),
      min_cells_per_sample = min(cell_counts$n_cells),
      genes_ranked_n = length(rank_stats),
      tested_KEGG_sets_n = length(tested_sets),
      direction = dplyr::if_else(NES > 0, "relapse_up", "relapse_down"),
      .before = 1
    )

  list(
    result = result,
    logcpm = if (identical(spec$candidate_id, "A05")) edgeR::cpm(y, log = TRUE, prior.count = 2) else NULL,
    cell_counts = cell_counts
  )
}

results <- vector("list", nrow(specs))
a05_logcpm <- NULL
for (i in seq_len(nrow(specs))) {
  cat("Running", specs$candidate_id[i], "-", specs$candidate_label[i], "\n")
  run <- run_candidate(specs[i, , drop = FALSE], i)
  results[[i]] <- run$result
  if (identical(specs$candidate_id[i], "A05")) a05_logcpm <- run$logcpm
}
all_results <- dplyr::bind_rows(results) |>
  dplyr::arrange(definition_changes, candidate_id, padj, dplyr::desc(abs(NES)))

target_summary <- all_results |>
  dplyr::filter(pathway_id %in% c("hsa03050", "hsa04137")) |>
  dplyr::mutate(
    global_FDR_lt_0.05 = padj < 0.05,
    global_FDR_lt_0.10 = padj < 0.10
  ) |>
  dplyr::arrange(definition_changes, candidate_id, pathway_id)

# Patient-direction audit for the recommended high-specificity candidate.
a05_targets <- target_summary |>
  dplyr::filter(candidate_id == "A05")
a05_pair_direction <- lapply(seq_len(nrow(a05_targets)), function(i) {
  genes <- strsplit(a05_targets$leadingEdge[[i]], ";", fixed = TRUE)[[1]]
  genes <- intersect(genes, rownames(a05_logcpm))
  score <- colMeans(a05_logcpm[genes, , drop = FALSE])
  tibble::tibble(
    pathway_id = a05_targets$pathway_id[[i]],
    pathway_name = a05_targets$pathway_name[[i]],
    GSEA_NES = a05_targets$NES[[i]],
    GSEA_nominal_P = a05_targets$pval[[i]],
    GSEA_global_FDR = a05_targets$padj[[i]],
    leading_edge_genes_n = length(genes),
    sample_id = names(score),
    leading_edge_mean_logCPM = as.numeric(score)
  ) |>
    dplyr::left_join(sample_info, by = "sample_id")
}) |>
  dplyr::bind_rows() |>
  dplyr::group_by(pathway_id, pathway_name, patient_pair) |>
  dplyr::summarise(
    diagnosis_score = leading_edge_mean_logCPM[condition_simple == "diagnosis"],
    relapse_score = leading_edge_mean_logCPM[condition_simple == "relapse"],
    paired_delta = relapse_score - diagnosis_score,
    leading_edge_genes_n = dplyr::first(leading_edge_genes_n),
    GSEA_NES = dplyr::first(GSEA_NES),
    GSEA_nominal_P = dplyr::first(GSEA_nominal_P),
    GSEA_global_FDR = dplyr::first(GSEA_global_FDR),
    .groups = "drop"
  ) |>
  dplyr::mutate(relapse_higher = paired_delta > 0)

write_csv_file(
  all_results,
  file.path(paths$tables, "publication_compatible_candidate_full_KEGG_all_results.csv")
)
write_csv_file(
  target_summary,
  file.path(paths$tables, "publication_compatible_candidate_full_KEGG_target_summary.csv")
)
write_csv_file(
  a05_pair_direction,
  file.path(paths$tables, "publication_compatible_A05_target_pair_direction.csv")
)

cat("Done. Candidates:", nrow(specs), "; full rows:", nrow(all_results), "\n")
print(target_summary |>
  dplyr::select(
    candidate_id, candidate_label, n_cells_total, min_cells_per_sample,
    pathway_name, NES, pval, padj, global_FDR_lt_0.05, global_FDR_lt_0.10
  ))
print(a05_pair_direction)
