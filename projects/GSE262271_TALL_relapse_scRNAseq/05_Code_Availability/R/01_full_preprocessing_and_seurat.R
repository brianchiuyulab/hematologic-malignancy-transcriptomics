options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grepl("^--file=", args)]
script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE) else normalizePath("gse262271_full_analysis.R", winslash = "/", mustWork = FALSE)
script_dir <- dirname(script_path)
trailing <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = "") {
  hit <- match(flag, trailing)
  if (!is.na(hit) && length(trailing) >= hit + 1) trailing[[hit + 1]] else default
}

default_project <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
project <- normalizePath(get_arg("--project", default_project), winslash = "/", mustWork = FALSE)
set.seed(20260704)

paths <- list(
  root = project,
  code = file.path(project, "code avalibility"),
  source = file.path(project, "source_metadata"),
  raw = file.path(project, "raw_data"),
  objects = file.path(project, "processed_objects"),
  results = file.path(project, "results"),
  figures = file.path(project, "figures"),
  logs = file.path(project, "logs")
)
invisible(lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))

log_file <- file.path(paths$logs, "pipeline.log")
log_message <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

load_or_stop <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(lapply(pkgs, library, character.only = TRUE))
}

load_or_stop(c("Seurat", "SeuratObject", "Matrix", "ggplot2", "dplyr", "tidyr", "readr", "patchwork", "harmony"))

theme_pub <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(color = "black"),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "grey70"),
      legend.key = ggplot2::element_blank()
    )
}

get_assay_matrix <- function(obj, assay = "RNA", layer = "data") {
  tryCatch(
    Seurat::GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) Seurat::GetAssayData(obj, assay = assay, slot = layer)
  )
}

source_notes <- c(
  "GSE262271 source paper:",
  "Kypraios et al. Identifying Candidate Gene Drivers Associated with Relapse in Pediatric T-Cell Acute Lymphoblastic Leukemia Using a Gene Co-Expression Network Approach. Cancers 2024. DOI: 10.3390/cancers16091667.",
  "The paper states that GSE262271 contains three paired pediatric T-ALL diagnosis-relapse scRNA-seq datasets.",
  "Supplementary information text maps samples as: Patient#1 M104 diagnosis / M127 relapse; Patient#2 M143 diagnosis with corticosteroid resistance / M148 relapse; Patient#3 M187 diagnosis / M187r relapse.",
  "The paper reports 10x Chromium Next GEM Single Cell 5' v2 with feature barcoding, Cell Ranger v6.0.0, Seurat v5.0.2, HTODemux default parameters, MT < 5%, and sample-specific lower nFeature cutoffs: M104/M127 >200, M143/M148 >300, M187/M187r >90.",
  "The paper states R scripts can be found at Peyronlab.github.io; in this run the public Peyronlab GitHub organization only exposed general repositories and no article-specific GSE262271 pipeline repository."
)
writeLines(source_notes, file.path(paths$source, "source_and_code_availability_notes.txt"))

sample_meta <- data.frame(
  library_id = c("lib_1", "lib_1", "lib_2", "lib_2", "lib_3", "lib_3"),
  geo_accession = c("GSM8162357", "GSM8162357", "GSM8162358", "GSM8162358", "GSM8162359", "GSM8162359"),
  file = c("GSM8162357_220415_lib_1.h5", "GSM8162357_220415_lib_1.h5",
           "GSM8162358_220415_lib_2.h5", "GSM8162358_220415_lib_2.h5",
           "GSM8162359_220415_lib_3.h5", "GSM8162359_220415_lib_3.h5"),
  patient_pair = c("P1", "P1", "P2", "P2", "P3", "P3"),
  hto_tag = c("M104", "M127", "M143", "M148", "M187", "M187r"),
  condition = c("diagnosis", "relapse", "diagnosis_corticoR", "relapse", "diagnosis", "relapse"),
  condition_simple = c("diagnosis", "relapse", "diagnosis", "relapse", "diagnosis", "relapse"),
  paper_label = c("M104 diag", "M127 rel", "M143 diag with cortico resistance", "M148 rel", "M187 diag", "M187r rel"),
  min_features_paper = c(200, 200, 300, 300, 90, 90)
)
readr::write_csv(sample_meta, file.path(paths$source, "sample_tag_metadata.csv"))

geo_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE262nnn/GSE262271/suppl/GSE262271_RAW.tar"
tar_file <- file.path(paths$raw, "GSE262271_RAW.tar")
if (!file.exists(tar_file)) {
  log_message("Downloading", geo_url)
  utils::download.file(geo_url, tar_file, mode = "wb", quiet = FALSE)
}
extract_dir <- file.path(paths$raw, "GSE262271_RAW_extracted")
dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
if (length(list.files(extract_dir, pattern = "\\.h5$", full.names = TRUE)) < 3) {
  log_message("Extracting", tar_file)
  utils::untar(tar_file, exdir = extract_dir)
}

h5_files <- file.path(extract_dir, unique(sample_meta$file))
if (!all(file.exists(h5_files))) stop("Missing H5 files: ", paste(h5_files[!file.exists(h5_files)], collapse = ", "), call. = FALSE)

raw_objs <- list()
demux_summaries <- list()
for (f in h5_files) {
  lib_file <- basename(f)
  lib_meta <- sample_meta[sample_meta$file == lib_file, ]
  lib_id <- lib_meta$library_id[1]
  log_message("Reading", lib_file)
  mats <- Seurat::Read10X_h5(f, use.names = TRUE, unique.features = TRUE)
  if (!is.list(mats) || !"Gene Expression" %in% names(mats) || !"Antibody Capture" %in% names(mats)) {
    stop("Expected Gene Expression and Antibody Capture matrices in ", lib_file, call. = FALSE)
  }
  obj <- Seurat::CreateSeuratObject(counts = mats[["Gene Expression"]], project = lib_id, min.cells = 3, min.features = 0)
  hto <- mats[["Antibody Capture"]][lib_meta$hto_tag, colnames(obj), drop = FALSE]
  obj[["HTO"]] <- Seurat::CreateAssayObject(counts = hto)
  obj$library_id <- lib_id
  obj$geo_accession <- lib_meta$geo_accession[1]
  obj$patient_pair <- lib_meta$patient_pair[1]
  obj <- Seurat::NormalizeData(obj, assay = "HTO", normalization.method = "CLR", margin = 2, verbose = FALSE)
  obj <- Seurat::HTODemux(obj, assay = "HTO", positive.quantile = 0.99, verbose = FALSE)
  obj$hto_tag <- as.character(obj$HTO_maxID)
  obj$sample_id <- sample_meta$hto_tag[match(obj$hto_tag, sample_meta$hto_tag)]
  obj$condition <- sample_meta$condition[match(obj$hto_tag, sample_meta$hto_tag)]
  obj$condition_simple <- sample_meta$condition_simple[match(obj$hto_tag, sample_meta$hto_tag)]
  obj$paper_label <- sample_meta$paper_label[match(obj$hto_tag, sample_meta$hto_tag)]
  obj$min_features_paper <- sample_meta$min_features_paper[match(obj$hto_tag, sample_meta$hto_tag)]
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-")

  demux_summaries[[lib_id]] <- as.data.frame(table(HTO_classification.global = obj$HTO_classification.global, HTO_maxID = obj$HTO_maxID)) |>
    dplyr::mutate(library_id = lib_id)

  png(file.path(paths$figures, paste0("HTO_demux_", lib_id, ".png")), width = 1600, height = 1000, res = 150)
  print(
    Seurat::RidgePlot(obj, assay = "HTO", features = rownames(obj[["HTO"]]), ncol = 2) /
      (ggplot2::ggplot(obj@meta.data, ggplot2::aes(x = HTO_classification.global, fill = HTO_maxID)) +
         ggplot2::geom_bar() +
         ggplot2::labs(title = paste("HTO demux", lib_id), x = NULL, y = "Cells") +
         theme_pub())
  )
  grDevices::dev.off()

  obj <- Seurat::RenameCells(obj, add.cell.id = lib_id)
  raw_objs[[lib_id]] <- obj
}

demux_summary <- dplyr::bind_rows(demux_summaries)
readr::write_csv(demux_summary, file.path(paths$results, "HTO_demux_summary.csv"))

qc_pre <- dplyr::bind_rows(lapply(raw_objs, function(o) {
  data.frame(
    cell = colnames(o),
    library_id = o$library_id,
    hto_classification = o$HTO_classification.global,
    hto_tag = o$hto_tag,
    sample_id = o$sample_id,
    patient_pair = o$patient_pair,
    condition = o$condition,
    condition_simple = o$condition_simple,
    nCount_RNA = o$nCount_RNA,
    nFeature_RNA = o$nFeature_RNA,
    percent.mt = o$percent.mt
  )
}))
readr::write_csv(qc_pre, file.path(paths$results, "qc_prefilter_cells.csv"))

qc_objs <- lapply(raw_objs, function(o) {
  keep <- o$HTO_classification.global == "Singlet" &
    !is.na(o$sample_id) &
    o$percent.mt < 5 &
    o$nFeature_RNA >= o$min_features_paper
  subset(o, cells = colnames(o)[keep])
})
qc_objs <- qc_objs[vapply(qc_objs, ncol, numeric(1)) > 0]
if (!length(qc_objs)) stop("No cells passed HTO singlet and QC filters.", call. = FALSE)

qc_post <- dplyr::bind_rows(lapply(qc_objs, function(o) {
  data.frame(
    cell = colnames(o),
    library_id = o$library_id,
    sample_id = o$sample_id,
    patient_pair = o$patient_pair,
    condition = o$condition,
    condition_simple = o$condition_simple,
    nCount_RNA = o$nCount_RNA,
    nFeature_RNA = o$nFeature_RNA,
    percent.mt = o$percent.mt
  )
}))
readr::write_csv(qc_post, file.path(paths$results, "qc_postfilter_cells.csv"))

qc_summary <- qc_pre |>
  dplyr::count(library_id, sample_id, patient_pair, condition_simple, name = "cells_prefilter") |>
  dplyr::full_join(
    qc_post |> dplyr::count(library_id, sample_id, patient_pair, condition_simple, name = "cells_postfilter"),
    by = c("library_id", "sample_id", "patient_pair", "condition_simple")
  ) |>
  dplyr::mutate(cells_postfilter = ifelse(is.na(cells_postfilter), 0, cells_postfilter))
readr::write_csv(qc_summary, file.path(paths$results, "qc_filter_summary.csv"))

png(file.path(paths$figures, "QC_violin_prefilter_postfilter.png"), width = 1800, height = 1200, res = 150)
pre_plot <- ggplot2::ggplot(qc_pre, ggplot2::aes(sample_id, nFeature_RNA, fill = condition_simple)) +
  ggplot2::geom_violin(scale = "width", trim = TRUE) +
  ggplot2::geom_hline(ggplot2::aes(yintercept = min_features_paper), data = sample_meta, linetype = 2) +
  ggplot2::labs(title = "Pre-filter nFeature_RNA", x = NULL, y = "nFeature_RNA") + theme_pub()
post_plot <- ggplot2::ggplot(qc_post, ggplot2::aes(sample_id, nFeature_RNA, fill = condition_simple)) +
  ggplot2::geom_violin(scale = "width", trim = TRUE) +
  ggplot2::labs(title = "Post-filter nFeature_RNA", x = NULL, y = "nFeature_RNA") + theme_pub()
mt_plot <- ggplot2::ggplot(qc_post, ggplot2::aes(sample_id, percent.mt, fill = condition_simple)) +
  ggplot2::geom_violin(scale = "width", trim = TRUE) +
  ggplot2::geom_hline(yintercept = 5, linetype = 2) +
  ggplot2::labs(title = "Post-filter percent.mt", x = NULL, y = "percent.mt") + theme_pub()
print(pre_plot / post_plot / mt_plot)
grDevices::dev.off()

obj <- if (length(qc_objs) == 1) qc_objs[[1]] else merge(qc_objs[[1]], y = qc_objs[-1])
if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
  obj <- SeuratObject::JoinLayers(obj, assay = "RNA")
}
DefaultAssay(obj) <- "RNA"

scale_factor <- median(obj$nCount_RNA)
log_message("Normalizing with LogNormalize scale.factor =", round(scale_factor, 2))
obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = scale_factor, verbose = FALSE)
obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
obj <- Seurat::ScaleData(obj, features = Seurat::VariableFeatures(obj), verbose = FALSE)
obj <- Seurat::RunPCA(obj, features = Seurat::VariableFeatures(obj), npcs = 50, verbose = FALSE)

png(file.path(paths$figures, "PCA_elbow_plot.png"), width = 1050, height = 750, res = 150)
print(Seurat::ElbowPlot(obj, ndims = 50))
grDevices::dev.off()
pca_sd <- obj[["pca"]]@stdev
pca_tbl <- data.frame(PC = seq_along(pca_sd), stdev = pca_sd, variance = pca_sd^2 / sum(pca_sd^2), cumulative_variance = cumsum(pca_sd^2 / sum(pca_sd^2)))
readr::write_csv(pca_tbl, file.path(paths$results, "pca_variance_table.csv"))

dims_use <- 1:30
obj <- harmony::RunHarmony(obj, group.by.vars = "library_id", reduction.use = "pca", dims.use = dims_use, theta = 2, verbose = FALSE)
obj <- Seurat::RunUMAP(obj, reduction = "harmony", dims = dims_use, reduction.name = "umap.harmony", verbose = FALSE)
obj <- Seurat::FindNeighbors(obj, reduction = "harmony", dims = dims_use, verbose = FALSE)
for (res in c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)) {
  obj <- Seurat::FindClusters(obj, resolution = res, verbose = FALSE)
}
obj$cluster_res0.6 <- as.character(obj$RNA_snn_res.0.6)
Idents(obj) <- "cluster_res0.6"

png(file.path(paths$figures, "UMAP_sample_condition_cluster.png"), width = 2100, height = 1500, res = 150)
print(
  Seurat::DimPlot(obj, reduction = "umap.harmony", group.by = "sample_id", raster = FALSE) +
    Seurat::DimPlot(obj, reduction = "umap.harmony", group.by = "condition_simple", raster = FALSE) +
    Seurat::DimPlot(obj, reduction = "umap.harmony", group.by = "patient_pair", raster = FALSE) +
    Seurat::DimPlot(obj, reduction = "umap.harmony", group.by = "cluster_res0.6", label = TRUE, repel = TRUE, raster = FALSE)
)
grDevices::dev.off()

saveRDS(obj, file.path(paths$objects, "GSE262271_seurat_qc_harmony_unannotated.rds"))
capture.output(sessionInfo(), file = file.path(paths$logs, paste0("sessionInfo_prepare_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")))
log_message("Prepared QC/Harmony unannotated Seurat object. Run gse262271_downstream_analysis.R for manual annotation, pathway summaries, and inferCNV inputs.")
quit(save = "no", status = 0)
