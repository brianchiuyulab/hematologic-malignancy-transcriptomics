options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = "") {
  hit <- match(flag, args)
  if (!is.na(hit) && length(args) >= hit + 1) args[[hit + 1]] else default
}

default_project <- file.path(
  "C:/Users/User/Desktop/Public dataset analysis/GEO dataset/2026/single cell",
  "BTZ \u76f8\u95dc",
  "GSE262271_TALL_relapse_scRNAseq"
)
project <- normalizePath(get_arg("--project", default_project), winslash = "/", mustWork = TRUE)

paths <- list(
  objects = file.path(project, "processed_objects"),
  results = file.path(project, "results"),
  figures = file.path(project, "figures"),
  logs = file.path(project, "logs")
)
invisible(lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))

user_lib <- Sys.getenv("R_LIBS_USER", unset = "C:/Users/User/Documents/R/win-library/4.4")
if (dir.exists(user_lib)) .libPaths(unique(c(normalizePath(user_lib, winslash = "/"), .libPaths())))

load_or_stop <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
  suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))
}
load_or_stop(c("Seurat", "SeuratObject", "ggplot2", "dplyr", "tibble", "patchwork"))

write_csv_file <- function(x, path) {
  utils::write.csv(as.data.frame(x), path, row.names = FALSE, fileEncoding = "UTF-8")
}

theme_umap <- function(base_size = 9) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(color = "black", size = base_size - 1),
      axis.title = ggplot2::element_text(color = "black"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = base_size + 1),
      legend.title = ggplot2::element_blank(),
      legend.key.height = grid::unit(0.35, "cm"),
      legend.key.width = grid::unit(0.35, "cm"),
      legend.text = ggplot2::element_text(size = base_size - 1)
    )
}

theme_dot <- function(base_size = 8) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
      axis.text.y = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0)
    )
}

join_layers_safe <- function(obj) {
  if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
    obj <- SeuratObject::JoinLayers(obj, assay = "RNA")
  }
  obj
}

obj_file <- file.path(paths$objects, "GSE262271_seurat_qc_harmony_annotated.rds")
if (!file.exists(obj_file)) stop("Missing annotated Seurat object: ", obj_file, call. = FALSE)
obj <- readRDS(obj_file)
obj <- join_layers_safe(obj)
Seurat::DefaultAssay(obj) <- "RNA"

required_cols <- c("sample_id", "patient_pair", "condition_simple", "cluster_res0.6", "major_annotation", "broad_annotation", "candidate_malignant_like")
missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
if (length(missing_cols)) stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)

major_map <- c(
  "B-cell-like_contaminant" = "B cell-like",
  "Cycling_TALL_blast_like" = "Cycling T-ALL-like",
  "Cytotoxic_T_or_NK_like" = "Cytotoxic T/NK-like",
  "ER_stress_TALL_like" = "ER-stress T-ALL-like",
  "Erythroid-like_contaminant" = "Erythroid-like",
  "Immature_HSC_like_TALL" = "Immature/HSC-like T-ALL",
  "Myeloid-like_contaminant" = "Myeloid-like",
  "pDC-like_contaminant" = "pDC-like",
  "Quiescent_CD44_high_TALL_like" = "CD44-high quiescent T-ALL-like",
  "TALL_blast_like" = "T-ALL blast-like",
  "T_lineage_like" = "T-lineage-like"
)
broad_map <- c(
  "B_cell_reference" = "B cell reference",
  "Cytotoxic_T_or_NK_reference_like" = "Cytotoxic T/NK reference-like",
  "Erythroid_reference" = "Erythroid reference",
  "Myeloid_reference" = "Myeloid reference",
  "pDC_reference" = "pDC reference",
  "Plasma_reference" = "Plasma reference",
  "T_lineage_candidate_malignant" = "T-lineage candidate malignant"
)

obj$major_annotation_publication <- unname(major_map[as.character(obj$major_annotation)])
obj$major_annotation_publication[is.na(obj$major_annotation_publication)] <- as.character(obj$major_annotation[is.na(obj$major_annotation_publication)])
obj$broad_annotation_publication <- unname(broad_map[as.character(obj$broad_annotation)])
obj$broad_annotation_publication[is.na(obj$broad_annotation_publication)] <- as.character(obj$broad_annotation[is.na(obj$broad_annotation_publication)])
obj$condition_publication <- ifelse(obj$condition_simple == "relapse", "Relapse", "Diagnosis")

emb <- as.data.frame(Seurat::Embeddings(obj, reduction = "umap.harmony"))
colnames(emb)[1:2] <- c("UMAP_1", "UMAP_2")
meta <- tibble::rownames_to_column(obj@meta.data, "cell")
plot_df <- cbind(meta, emb[meta$cell, c("UMAP_1", "UMAP_2")])

cnv_file <- file.path(paths$results, "infercnv_nonTref_cell_cnv_burden.csv")
plot_df$strict_cancer_like <- ifelse(plot_df$candidate_malignant_like, "Candidate T-lineage", "Reference-like / non-candidate")
if (file.exists(cnv_file)) {
  cnv <- utils::read.csv(cnv_file, stringsAsFactors = FALSE, check.names = FALSE)
  cell_col <- intersect(c("cell", "barcode", "cell_id", "cell_name"), colnames(cnv))
  if (!length(cell_col)) cell_col <- colnames(cnv)[[1]]
  score_col <- intersect(c("infercnv_mean_abs_log2", "mean_abs_log2"), colnames(cnv))
  if (!length(score_col)) score_col <- grep("mean.*abs.*log2|abs.*log2", colnames(cnv), ignore.case = TRUE, value = TRUE)
  if (length(score_col)) {
    cnv_score <- cnv[[score_col[[1]]]]
    names(cnv_score) <- cnv[[cell_col[[1]]]]
    plot_df$infercnv_mean_abs_log2 <- as.numeric(cnv_score[plot_df$cell])
    cand_scores <- plot_df$infercnv_mean_abs_log2[plot_df$candidate_malignant_like & is.finite(plot_df$infercnv_mean_abs_log2)]
    if (length(cand_scores) > 100) {
      q50 <- stats::quantile(cand_scores, 0.50, na.rm = TRUE)
      q75 <- stats::quantile(cand_scores, 0.75, na.rm = TRUE)
      plot_df$strict_cancer_like[plot_df$candidate_malignant_like & plot_df$infercnv_mean_abs_log2 >= q50] <- "CNV-high candidate top 50%"
      plot_df$strict_cancer_like[plot_df$candidate_malignant_like & plot_df$infercnv_mean_abs_log2 >= q75] <- "CNV-high candidate top 25%"
    }
  }
}
plot_df$strict_cancer_like <- factor(
  plot_df$strict_cancer_like,
  levels = c("Reference-like / non-candidate", "Candidate T-lineage", "CNV-high candidate top 50%", "CNV-high candidate top 25%")
)

major_cols <- c(
  "B cell-like" = "#4E79A7",
  "Cycling T-ALL-like" = "#D95F02",
  "Cytotoxic T/NK-like" = "#59A14F",
  "ER-stress T-ALL-like" = "#B07AA1",
  "Erythroid-like" = "#E15759",
  "Immature/HSC-like T-ALL" = "#00A0B0",
  "Myeloid-like" = "#9C755F",
  "pDC-like" = "#EDC948",
  "CD44-high quiescent T-ALL-like" = "#7B5EA7",
  "T-ALL blast-like" = "#F28E2B",
  "T-lineage-like" = "#76B7B2"
)
broad_cols <- c(
  "T-lineage candidate malignant" = "#C7443E",
  "B cell reference" = "#4E79A7",
  "Cytotoxic T/NK reference-like" = "#59A14F",
  "Erythroid reference" = "#E15759",
  "Myeloid reference" = "#9C755F",
  "pDC reference" = "#EDC948",
  "Plasma reference" = "#B07AA1"
)
strict_cols <- c(
  "Reference-like / non-candidate" = "grey80",
  "Candidate T-lineage" = "#D6A33A",
  "CNV-high candidate top 50%" = "#C7443E",
  "CNV-high candidate top 25%" = "#7A1F1D"
)
condition_cols <- c("Diagnosis" = "#3C6E71", "Relapse" = "#B2433F")

base_scatter <- function(df, color_col, colors, title) {
  ggplot2::ggplot(df, ggplot2::aes(x = UMAP_1, y = UMAP_2, color = .data[[color_col]])) +
    ggplot2::geom_point(size = 0.12, alpha = 0.78, stroke = 0) +
    ggplot2::coord_equal() +
    ggplot2::scale_color_manual(values = colors, drop = FALSE) +
    ggplot2::labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 2.4, alpha = 1), ncol = 1)) +
    theme_umap(9)
}

p_major <- base_scatter(plot_df, "major_annotation_publication", major_cols, "Marker-based major annotation")
p_broad <- base_scatter(plot_df, "broad_annotation_publication", broad_cols, "Broad malignant/reference grouping")
p_strict <- base_scatter(plot_df, "strict_cancer_like", strict_cols, "CNV-strict candidate malignant-like cells")
p_condition <- base_scatter(plot_df, "condition_publication", condition_cols, "Diagnosis vs relapse")

umap_combined <- (p_major + p_broad) / (p_strict + p_condition) +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "right")
ggplot2::ggsave(file.path(paths$figures, "UMAP_annotation_publication_QC.png"), umap_combined, width = 12.5, height = 8.5, dpi = 300)
ggplot2::ggsave(file.path(paths$figures, "UMAP_strict_CNV_candidate_malignant_like.png"), p_strict + ggplot2::theme(legend.position = "right"), width = 6.8, height = 5.0, dpi = 300)

cluster_centers <- plot_df |>
  dplyr::group_by(cluster_res0.6) |>
  dplyr::summarise(UMAP_1 = stats::median(UMAP_1), UMAP_2 = stats::median(UMAP_2), .groups = "drop")
p_cluster <- ggplot2::ggplot(plot_df, ggplot2::aes(UMAP_1, UMAP_2, color = cluster_res0.6)) +
  ggplot2::geom_point(size = 0.10, alpha = 0.75, stroke = 0) +
  ggplot2::geom_label(
    data = cluster_centers,
    ggplot2::aes(x = UMAP_1, y = UMAP_2, label = cluster_res0.6),
    inherit.aes = FALSE,
    size = 2.8,
    linewidth = 0.15,
    fill = "white",
    alpha = 0.88
  ) +
  ggplot2::coord_equal() +
  ggplot2::labs(title = "Seurat SNN clusters at resolution 0.6", x = "UMAP 1", y = "UMAP 2") +
  theme_umap(9) +
  ggplot2::theme(legend.position = "none")
ggplot2::ggsave(file.path(paths$figures, "UMAP_cluster_res0.6_publication_QC.png"), p_cluster, width = 6.8, height = 5.0, dpi = 300)

lineage_features <- c(
  "CD3D", "CD3E", "TRAC", "CD7", "CD2",
  "NKG7", "GNLY", "PRF1", "KLRD1",
  "MS4A1", "CD79A", "CD74",
  "LYZ", "S100A8", "S100A9", "FCGR3A",
  "LILRA4", "IL3RA", "GZMB", "IRF7",
  "HBB", "HBA1", "HBA2", "ALAS2", "GYPA"
)
tall_state_features <- c(
  "CD7", "CD3D", "TRAC", "PTCRA", "DNTT", "HES1", "NOTCH1",
  "CD34", "KIT", "SPINK2", "MEIS1", "PROM1",
  "CD44", "IL7R", "BCL2",
  "MKI67", "TOP2A", "TYMS", "STMN1", "HMGB2",
  "XBP1", "HSPA5", "DDIT3", "ATF4"
)
lineage_features <- intersect(lineage_features, rownames(obj))
tall_state_features <- intersect(tall_state_features, rownames(obj))

if (length(lineage_features) >= 5) {
  p_dot_lineage <- Seurat::DotPlot(
    obj,
    features = lineage_features,
    group.by = "major_annotation_publication",
    dot.scale = 5
  ) +
    ggplot2::scale_color_gradient2(low = "grey88", mid = "#8EA8C3", high = "#12355B", midpoint = 0) +
    ggplot2::labs(title = "Lineage marker check") +
    theme_dot(8)
  ggplot2::ggsave(file.path(paths$figures, "annotation_lineage_marker_dotplot_publication.png"), p_dot_lineage, width = 11.5, height = 5.2, dpi = 300)
}

if (length(tall_state_features) >= 5) {
  p_dot_state <- Seurat::DotPlot(
    obj,
    features = tall_state_features,
    group.by = "major_annotation_publication",
    dot.scale = 5
  ) +
    ggplot2::scale_color_gradient2(low = "grey88", mid = "#B8A1D9", high = "#5B2A86", midpoint = 0) +
    ggplot2::labs(title = "T-ALL state marker check") +
    theme_dot(8)
  ggplot2::ggsave(file.path(paths$figures, "annotation_TALL_state_marker_dotplot_publication.png"), p_dot_state, width = 11.5, height = 5.2, dpi = 300)
}

cluster_annotation <- file.path(paths$results, "cluster_annotation_res0.6.csv")
top_marker_file <- file.path(paths$results, "cluster_top20_markers_res0.6_fast_review.csv")
if (file.exists(cluster_annotation) && file.exists(top_marker_file)) {
  ann <- utils::read.csv(cluster_annotation, stringsAsFactors = FALSE, check.names = FALSE)
  top <- utils::read.csv(top_marker_file, stringsAsFactors = FALSE, check.names = FALSE)
  gene_col <- if ("gene" %in% colnames(top)) "gene" else intersect(c("features", "feature"), colnames(top))[1]
  if (!is.na(gene_col)) {
    top5 <- top |>
      dplyr::group_by(cluster) |>
      dplyr::slice_head(n = 5) |>
      dplyr::summarise(top5_markers = paste(.data[[gene_col]], collapse = ";"), .groups = "drop")
    ann$cluster <- ann$cluster_res0.6
    out <- dplyr::left_join(ann, top5, by = "cluster")
    write_csv_file(out, file.path(paths$results, "annotation_cluster_top5_marker_audit.csv"))
  }
}

annotation_counts <- plot_df |>
  dplyr::count(major_annotation_publication, broad_annotation_publication, condition_publication, sample_id, name = "n_cells") |>
  dplyr::arrange(major_annotation_publication, sample_id)
write_csv_file(annotation_counts, file.path(paths$results, "annotation_publication_cell_counts_by_sample.csv"))

report <- c(
  "# Annotation QC publication figures",
  "",
  "- Annotation was generated from Seurat SNN/KNN graph clustering after Harmony correction, not SingleR.",
  "- `cluster_res0.6` is the active clustering resolution used for manual marker-panel annotation.",
  "- Candidate malignant-like cells are the broad T-lineage candidate compartment; stricter sensitivity uses non-T-reference inferCNV burden top 50% and top 25% within that compartment.",
  "- New publication-style figures avoid overlaid UMAP text labels and use focused marker panels instead of a very wide all-marker DotPlot.",
  "",
  "## Files",
  "- figures/UMAP_annotation_publication_QC.png",
  "- figures/UMAP_strict_CNV_candidate_malignant_like.png",
  "- figures/UMAP_cluster_res0.6_publication_QC.png",
  "- figures/annotation_lineage_marker_dotplot_publication.png",
  "- figures/annotation_TALL_state_marker_dotplot_publication.png",
  "- results/annotation_cluster_top5_marker_audit.csv",
  "- results/annotation_publication_cell_counts_by_sample.csv"
)
writeLines(report, file.path(paths$results, "ANNOTATION_QC_PUBLICATION_REPORT.md"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(paths$logs, paste0("sessionInfo_annotation_qc_publication_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")))
