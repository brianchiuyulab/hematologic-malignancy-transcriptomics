options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = "") {
  hit <- match(flag, args)
  if (!is.na(hit) && length(args) >= hit + 1) args[[hit + 1]] else default
}

project <- normalizePath(
  get_arg("--project", "C:/Users/User/Desktop/Public dataset analysis/GEO dataset/2026/single cell/BTZ \u76f8\u95dc/GSE262271_TALL_relapse_scRNAseq"),
  winslash = "/",
  mustWork = TRUE
)
paths <- list(
  objects = file.path(project, "processed_objects"),
  results = file.path(project, "results"),
  figures = file.path(project, "figures"),
  logs = file.path(project, "logs"),
  infercnv = file.path(project, "infercnv_inputs")
)
invisible(lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))

log_message <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = file.path(paths$logs, "downstream.log"), append = TRUE)
}

load_or_stop <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(lapply(pkgs, library, character.only = TRUE))
}
load_or_stop(c("Seurat", "SeuratObject", "Matrix", "ggplot2", "dplyr", "tidyr", "readr", "patchwork", "tibble"))

theme_pub <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(color = "black"),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "grey70"),
      legend.key = ggplot2::element_blank()
    )
}

join_layers_safe <- function(obj) {
  if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
    obj <- SeuratObject::JoinLayers(obj, assay = "RNA")
  }
  obj
}

get_assay_matrix <- function(obj, assay = "RNA", layer = "data") {
  tryCatch(
    Seurat::GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) {
      obj <<- join_layers_safe(obj)
      Seurat::GetAssayData(obj, assay = assay, layer = layer)
    }
  )
}

obj_file <- file.path(paths$objects, "GSE262271_seurat_qc_harmony_unannotated.rds")
if (!file.exists(obj_file)) stop("Missing unannotated Seurat object: ", obj_file, call. = FALSE)
log_message("Reading", obj_file)
obj <- readRDS(obj_file)
obj <- join_layers_safe(obj)
DefaultAssay(obj) <- "RNA"
if (!"cluster_res0.6" %in% colnames(obj@meta.data)) {
  obj$cluster_res0.6 <- as.character(obj$RNA_snn_res.0.6)
}
Idents(obj) <- "cluster_res0.6"

marker_panels <- list(
  T_lineage = c("CD3D", "CD3E", "TRAC", "CD2", "CD5", "CD7", "LCK"),
  TALL_blast_signature = c("CD7", "CD99", "HES4", "SOX4", "STMN1", "ARMH1", "DNTT", "PTCRA", "TCF7", "LEF1", "IL7R"),
  Immature_HSC_like = c("CD34", "PROM1", "KIT", "MEIS1", "HLF", "SPINK2"),
  Proliferation = c("MKI67", "TOP2A", "UBE2C", "HMGB2", "PCNA", "TYMS"),
  Quiescent_CD44 = c("CD44", "CXCR4", "BCL2", "KLF2", "KLF6", "JUN", "FOS"),
  Cytotoxic_T_NK = c("NKG7", "GNLY", "PRF1", "GZMB", "KLRD1", "FCGR3A"),
  Myeloid = c("LYZ", "LST1", "S100A8", "S100A9", "CTSS", "FCN1"),
  B_cell = c("MS4A1", "CD79A", "CD79B", "CD74", "BANK1"),
  Plasma = c("MZB1", "XBP1", "JCHAIN", "IGHG1", "IGKC"),
  pDC = c("CLEC4C", "LILRA4", "IL3RA", "TCF4", "IRF7", "IRF8", "PLD4", "SCT"),
  ER_Stress = c("HSPA5", "HYOU1", "DNAJB9", "DDIT3", "DERL3", "TXNDC5", "PDIA4", "SDF2L1", "HERPUD1", "ATF4", "XBP1"),
  Erythroid = c("HBB", "HBA1", "HBA2", "ALAS2", "AHSP")
)
marker_panels <- lapply(marker_panels, function(g) intersect(unique(g), rownames(obj)))
marker_panels <- marker_panels[lengths(marker_panels) > 0]
readr::write_csv(
  data.frame(panel = names(marker_panels), genes_present = vapply(marker_panels, paste, collapse = ";", FUN.VALUE = character(1))),
  file.path(paths$results, "annotation_marker_panels_used.csv")
)

norm_data <- get_assay_matrix(obj, "RNA", "data")
for (nm in names(marker_panels)) {
  obj[[paste0(nm, "_score")]] <- Matrix::colMeans(norm_data[marker_panels[[nm]], , drop = FALSE])
}

markers_file <- file.path(paths$results, "cluster_markers_res0.6_fast_review.csv")
if (file.exists(markers_file)) {
  log_message("Reusing existing fast marker table", markers_file)
  all_markers <- readr::read_csv(markers_file, show_col_types = FALSE)
} else {
  log_message("Running fast Seurat FindAllMarkers for annotation review")
  all_markers <- Seurat::FindAllMarkers(
    obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.5,
    test.use = "wilcox",
    max.cells.per.ident = 800,
    random.seed = 20260704,
    verbose = FALSE
  )
  readr::write_csv(all_markers, markers_file)
}
top_markers <- all_markers |>
  dplyr::group_by(cluster) |>
  dplyr::slice_max(order_by = avg_log2FC, n = 20, with_ties = FALSE) |>
  dplyr::ungroup()
readr::write_csv(top_markers, file.path(paths$results, "cluster_top20_markers_res0.6_fast_review.csv"))

score_cols <- paste0(names(marker_panels), "_score")
cluster_scores <- obj@meta.data |>
  dplyr::group_by(cluster_res0.6) |>
  dplyr::summarise(dplyr::across(dplyr::all_of(score_cols), mean), n_cells = dplyr::n(), .groups = "drop")
readr::write_csv(cluster_scores, file.path(paths$results, "cluster_marker_scores_res0.6.csv"))

annotate_cluster <- function(row) {
  vals <- as.numeric(row[score_cols])
  names(vals) <- sub("_score$", "", score_cols)
  non_t_names <- intersect(c("Myeloid", "B_cell", "Plasma", "pDC", "Erythroid", "Cytotoxic_T_NK"), names(vals))
  max_non_t <- max(vals[non_t_names], na.rm = TRUE)
  t_signal <- max(vals[c("TALL_blast_signature", "T_lineage")], na.rm = TRUE)
  if (!is.na(vals["Erythroid"]) && vals["Erythroid"] > t_signal && vals["Erythroid"] > 0.15) return("Erythroid-like_contaminant")
  if (!is.na(vals["Myeloid"]) && vals["Myeloid"] > t_signal && vals["Myeloid"] > 0.15) return("Myeloid-like_contaminant")
  if (!is.na(vals["pDC"]) && vals["pDC"] > t_signal && vals["pDC"] > 0.35) return("pDC-like_contaminant")
  if (!is.na(vals["B_cell"]) && vals["B_cell"] > t_signal && vals["B_cell"] > 0.15) return("B-cell-like_contaminant")
  if (!is.na(vals["ER_Stress"]) && vals["ER_Stress"] > 0.80 && vals["TALL_blast_signature"] > 0.50) return("ER_stress_TALL_like")
  if (!is.na(vals["Plasma"]) && vals["Plasma"] > (t_signal + 0.30) && vals["Plasma"] > 0.50) return("Plasma-like_contaminant")
  if (!is.na(vals["Cytotoxic_T_NK"]) && vals["Cytotoxic_T_NK"] > vals["TALL_blast_signature"] && vals["Cytotoxic_T_NK"] > 0.2) return("Cytotoxic_T_or_NK_like")
  if (!is.na(vals["Proliferation"]) && vals["Proliferation"] > stats::median(cluster_scores$Proliferation_score, na.rm = TRUE) && vals["TALL_blast_signature"] > 0.03) return("Cycling_TALL_blast_like")
  if (!is.na(vals["Quiescent_CD44"]) && vals["Quiescent_CD44"] > stats::median(cluster_scores$Quiescent_CD44_score, na.rm = TRUE) && vals["TALL_blast_signature"] > 0.03) return("Quiescent_CD44_high_TALL_like")
  if (!is.na(vals["Immature_HSC_like"]) && vals["Immature_HSC_like"] > 0.03 && vals["TALL_blast_signature"] > 0.02) return("Immature_HSC_like_TALL")
  if (!is.na(vals["TALL_blast_signature"]) && vals["TALL_blast_signature"] >= 0.02 && t_signal >= max_non_t) return("TALL_blast_like")
  if (!is.na(vals["T_lineage"]) && vals["T_lineage"] >= 0.02 && t_signal >= max_non_t) return("T_lineage_like")
  paste0("Unassigned_top_", names(vals)[which.max(vals)])
}

cluster_annotation <- cluster_scores
cluster_annotation$major_annotation <- apply(cluster_scores, 1, annotate_cluster)
cluster_annotation$broad_annotation <- dplyr::case_when(
  cluster_annotation$major_annotation == "B-cell-like_contaminant" ~ "B_cell_reference",
  cluster_annotation$major_annotation == "Myeloid-like_contaminant" ~ "Myeloid_reference",
  cluster_annotation$major_annotation == "Erythroid-like_contaminant" ~ "Erythroid_reference",
  cluster_annotation$major_annotation == "pDC-like_contaminant" ~ "pDC_reference",
  cluster_annotation$major_annotation == "Plasma-like_contaminant" ~ "Plasma_reference",
  cluster_annotation$major_annotation == "Cytotoxic_T_or_NK_like" ~ "Cytotoxic_T_or_NK_reference_like",
  TRUE ~ "T_lineage_candidate_malignant"
)
readr::write_csv(cluster_annotation, file.path(paths$results, "cluster_annotation_res0.6.csv"))
map <- setNames(cluster_annotation$major_annotation, cluster_annotation$cluster_res0.6)
obj$major_annotation <- unname(map[obj$cluster_res0.6])
obj$major_annotation[is.na(obj$major_annotation)] <- "Unassigned"
broad_map <- setNames(cluster_annotation$broad_annotation, cluster_annotation$cluster_res0.6)
obj$broad_annotation <- unname(broad_map[obj$cluster_res0.6])
obj$broad_annotation[is.na(obj$broad_annotation)] <- "Unassigned"
obj$candidate_malignant_like <- obj$broad_annotation == "T_lineage_candidate_malignant"
obj$infercnv_group <- dplyr::case_when(
  obj$broad_annotation == "T_lineage_candidate_malignant" ~ "obs_T_lineage_candidate_malignant",
  obj$broad_annotation == "B_cell_reference" ~ "ref_B_cell",
  obj$broad_annotation == "Myeloid_reference" ~ "ref_Myeloid",
  obj$broad_annotation == "Erythroid_reference" ~ "ref_Erythroid",
  obj$broad_annotation == "pDC_reference" ~ "ref_pDC",
  obj$broad_annotation == "Plasma_reference" ~ "ref_Plasma",
  obj$broad_annotation == "Cytotoxic_T_or_NK_reference_like" ~ "ref_Cytotoxic_T_or_NK",
  TRUE ~ "exclude_Unassigned"
)

t_lineage_clusters <- cluster_annotation |>
  dplyr::filter(broad_annotation == "T_lineage_candidate_malignant") |>
  dplyr::select(cluster_res0.6, n_cells, major_annotation, broad_annotation)
readr::write_csv(t_lineage_clusters, file.path(paths$results, "t_lineage_candidate_clusters_res0.6.csv"))
readr::write_csv(as.data.frame(table(obj$broad_annotation, obj$sample_id), stringsAsFactors = FALSE),
                 file.path(paths$results, "broad_annotation_by_sample_cell_counts.csv"))

features_for_dotplot <- unique(unlist(marker_panels))
features_for_dotplot <- features_for_dotplot[features_for_dotplot %in% rownames(obj)]
png(file.path(paths$figures, "annotation_marker_dotplot_by_cluster_res0.6.png"), width = 2400, height = 1500, res = 150)
print(Seurat::DotPlot(obj, features = features_for_dotplot, group.by = "cluster_res0.6") + Seurat::RotatedAxis() + theme_pub(8))
grDevices::dev.off()
png(file.path(paths$figures, "annotation_marker_dotplot_by_broad_annotation.png"), width = 2400, height = 1200, res = 150)
print(Seurat::DotPlot(obj, features = features_for_dotplot, group.by = "broad_annotation") + Seurat::RotatedAxis() + theme_pub(8))
grDevices::dev.off()
png(file.path(paths$figures, "annotation_marker_dotplot_by_major_annotation.png"), width = 2400, height = 1500, res = 150)
print(Seurat::DotPlot(obj, features = features_for_dotplot, group.by = "major_annotation") + Seurat::RotatedAxis() + theme_pub(8))
grDevices::dev.off()

png(file.path(paths$figures, "UMAP_annotation_malignant_like.png"), width = 1800, height = 750, res = 150)
print(
  Seurat::DimPlot(obj, reduction = "umap.harmony", group.by = "major_annotation", label = TRUE, repel = TRUE, raster = FALSE) +
    Seurat::DimPlot(obj, reduction = "umap.harmony", group.by = "broad_annotation", label = TRUE, repel = TRUE, raster = FALSE) +
    Seurat::DimPlot(obj, reduction = "umap.harmony", group.by = "candidate_malignant_like", raster = FALSE)
)
grDevices::dev.off()

proteasome_core <- c("PSMA1","PSMA2","PSMA3","PSMA4","PSMA5","PSMA6","PSMA7",
                     "PSMB1","PSMB2","PSMB3","PSMB4","PSMB5","PSMB6","PSMB7","PSMB8","PSMB9","PSMB10",
                     "PSMC1","PSMC2","PSMC3","PSMC4","PSMC5","PSMC6",
                     "PSMD1","PSMD2","PSMD3","PSMD4","PSMD5","PSMD6","PSMD7","PSMD8","PSMD9","PSMD10","PSMD11","PSMD12","PSMD13","PSMD14",
                     "PSME1","PSME2","PSME3")
mitophagy_core <- c("PINK1","PRKN","PARK7","BNIP3","BNIP3L","FUNDC1","BCL2L13",
                    "OPTN","SQSTM1","CALCOCO2","TBK1","TAX1BP1","NBR1",
                    "MAP1LC3A","MAP1LC3B","MAP1LC3C","GABARAPL1","GABARAPL2",
                    "ATG5","ATG7","ATG12","ATG13","BECN1","ULK1","ULK2","MFN1","MFN2","DNM1L")
pathway_sets <- list(
  proteasome_core = intersect(proteasome_core, rownames(obj)),
  mitophagy_core = intersect(mitophagy_core, rownames(obj))
)
readr::write_csv(data.frame(pathway = names(pathway_sets), genes = vapply(pathway_sets, paste, collapse = ";", FUN.VALUE = character(1))),
                 file.path(paths$results, "pathway_gene_sets_used.csv"))
for (nm in names(pathway_sets)) {
  obj[[paste0(nm, "_avg_expr")]] <- Matrix::colMeans(norm_data[pathway_sets[[nm]], , drop = FALSE])
}

if ("STAT5B" %in% rownames(obj)) {
  obj$STAT5B_log_norm <- as.numeric(norm_data["STAT5B", ])
  obj$STAT5B_detected <- obj$STAT5B_log_norm > 0
} else {
  obj$STAT5B_log_norm <- NA_real_
  obj$STAT5B_detected <- NA
}

candidate_cells <- rownames(obj@meta.data)[obj$candidate_malignant_like]
candidate_meta <- obj@meta.data[candidate_cells, c("library_id","sample_id","patient_pair","condition","condition_simple","cluster_res0.6","major_annotation","broad_annotation","infercnv_group","nCount_RNA","nFeature_RNA","percent.mt","STAT5B_log_norm","STAT5B_detected","proteasome_core_avg_expr","mitophagy_core_avg_expr")]
readr::write_csv(tibble::rownames_to_column(candidate_meta, "cell"), file.path(paths$results, "candidate_malignant_like_cells_metadata.csv"))

counts <- get_assay_matrix(obj, "RNA", "counts")
infercnv_cells <- rownames(obj@meta.data)[obj$infercnv_group != "exclude_Unassigned"]
infercnv_annotations <- data.frame(cell = infercnv_cells, group = obj$infercnv_group[infercnv_cells])
readr::write_tsv(infercnv_annotations, file.path(paths$infercnv, "infercnv_cell_annotations.txt"), col_names = FALSE)
readr::write_csv(as.data.frame(table(infercnv_annotations$group), stringsAsFactors = FALSE),
                 file.path(paths$infercnv, "infercnv_group_cell_counts.csv"))

infercnv_counts <- counts[, infercnv_cells, drop = FALSE]
saveRDS(infercnv_counts, file.path(paths$infercnv, "infercnv_raw_counts_matrix.rds"))
if (requireNamespace("AnnotationDbi", quietly = TRUE) && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  genes <- rownames(infercnv_counts)
  suppressMessages({
    chr <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = genes, column = "CHR", keytype = "SYMBOL", multiVals = "first")
    start <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = genes, column = "CHRLOC", keytype = "SYMBOL", multiVals = "first")
    end <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = genes, column = "CHRLOCEND", keytype = "SYMBOL", multiVals = "first")
  })
  gene_order <- data.frame(
    gene = genes,
    chr = paste0("chr", chr),
    start = suppressWarnings(as.numeric(start)),
    end = suppressWarnings(as.numeric(end)),
    stringsAsFactors = FALSE
  )
  gene_order <- gene_order[!is.na(gene_order$chr) & gene_order$chr != "chrNA" &
                             !is.na(gene_order$start) & !is.na(gene_order$end), ]
  gene_order$start <- abs(gene_order$start)
  gene_order$end <- abs(gene_order$end)
  gene_order <- gene_order |>
    dplyr::mutate(start2 = pmin(start, end), end2 = pmax(start, end), start = start2, end = end2) |>
    dplyr::select(gene, chr, start, end)
  gene_order <- gene_order[gene_order$chr %in% paste0("chr", c(1:22, "X", "Y")), ]
  readr::write_tsv(gene_order, file.path(paths$infercnv, "infercnv_gene_order_from_orgHsEgDb.txt"), col_names = FALSE)
} else {
  writeLines("AnnotationDbi/org.Hs.eg.db not available; create gene_order_file before running inferCNV.",
             file.path(paths$infercnv, "gene_order_not_created.txt"))
}

infercnv_runner <- c(
  "options(stringsAsFactors = FALSE)",
  "user_lib <- Sys.getenv('R_LIBS_USER', unset = 'C:/Users/User/Documents/R/win-library/4.4')",
  "if (dir.exists(user_lib)) .libPaths(unique(c(normalizePath(user_lib, winslash = '/'), .libPaths())))",
  "configure_jags <- function() {",
  "  candidates <- unique(c(Sys.getenv('JAGS_HOME', unset = NA_character_), 'C:/Users/User/Documents/JAGS/JAGS-4.3.2_extracted', 'C:/Program Files/JAGS/JAGS-4.3.2', 'C:/Program Files (x86)/JAGS/JAGS-4.3.2'))",
  "  candidates <- candidates[!is.na(candidates) & nzchar(candidates) & dir.exists(candidates)]",
  "  if (!length(candidates)) return(invisible(FALSE))",
  "  jags_home <- normalizePath(candidates[[1]], winslash = '/', mustWork = TRUE)",
  "  Sys.setenv(JAGS_HOME = jags_home)",
  "  jags_paths <- file.path(jags_home, c('x64/bin', 'bin', 'x64/modules', 'modules'))",
  "  jags_paths <- jags_paths[dir.exists(jags_paths)]",
  "  if (length(jags_paths)) {",
  "    current_path <- strsplit(Sys.getenv('PATH'), .Platform$path.sep, fixed = TRUE)[[1]]",
  "    Sys.setenv(PATH = paste(unique(c(jags_paths, current_path)), collapse = .Platform$path.sep))",
  "  }",
  "  invisible(TRUE)",
  "}",
  "configure_jags()",
  "input_dir <- normalizePath(dirname(sub('^--file=', '', commandArgs(trailingOnly = FALSE)[grep('^--file=', commandArgs(trailingOnly = FALSE))[1]])), winslash = '/', mustWork = TRUE)",
  "if (!requireNamespace('infercnv', quietly = TRUE)) stop('Package infercnv is not installed. Install infercnv, then rerun this script.', call. = FALSE)",
  "counts <- readRDS(file.path(input_dir, 'infercnv_raw_counts_matrix.rds'))",
  "annot_file <- file.path(input_dir, 'infercnv_cell_annotations.txt')",
  "gene_order_file <- file.path(input_dir, 'infercnv_gene_order_from_orgHsEgDb.txt')",
  "annotations <- read.delim(annot_file, header = FALSE, stringsAsFactors = FALSE)",
  "ref_groups <- unique(annotations$V2[grepl('^ref_', annotations$V2)])",
  "out_dir <- file.path(input_dir, 'infercnv_output')",
  "dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)",
  "infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix = counts, annotations_file = annot_file, delim = '\\t', gene_order_file = gene_order_file, ref_group_names = ref_groups)",
  "infercnv_obj <- infercnv::run(infercnv_obj, cutoff = 0.1, out_dir = out_dir, cluster_by_groups = TRUE, denoise = TRUE, HMM = FALSE, num_threads = max(1, parallel::detectCores() - 1))",
  "saveRDS(infercnv_obj, file.path(out_dir, 'infercnv_obj.rds'))"
)
writeLines(infercnv_runner, file.path(paths$infercnv, "run_infercnv_if_installed.R"))

sample_summary <- candidate_meta |>
  dplyr::group_by(patient_pair, sample_id, condition_simple, condition) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    mean_STAT5B = mean(STAT5B_log_norm, na.rm = TRUE),
    median_STAT5B = median(STAT5B_log_norm, na.rm = TRUE),
    pct_STAT5B_detected = mean(STAT5B_detected, na.rm = TRUE) * 100,
    mean_proteasome_core = mean(proteasome_core_avg_expr, na.rm = TRUE),
    mean_mitophagy_core = mean(mitophagy_core_avg_expr, na.rm = TRUE),
    .groups = "drop"
)

pseudo_rows <- lapply(split(candidate_cells, obj$sample_id[candidate_cells]), function(cells) {
  sample_id <- unique(obj$sample_id[cells])
  total_counts <- sum(Matrix::colSums(counts[, cells, drop = FALSE]))
  stat5b_count <- if ("STAT5B" %in% rownames(counts)) sum(counts["STAT5B", cells]) else NA_real_
  data.frame(sample_id = sample_id, STAT5B_pseudobulk_logCPM = log1p(stat5b_count / total_counts * 1e6), total_counts = total_counts)
})
pseudo_df <- dplyr::bind_rows(pseudo_rows)
sample_summary <- dplyr::left_join(sample_summary, pseudo_df, by = "sample_id")
readr::write_csv(sample_summary, file.path(paths$results, "sample_level_candidate_malignant_like_summary.csv"))

paired_compare <- function(df, metric) {
  wide <- df |>
    dplyr::select(patient_pair, condition_simple, dplyr::all_of(metric)) |>
    tidyr::pivot_wider(names_from = condition_simple, values_from = dplyr::all_of(metric)) |>
    dplyr::mutate(delta_relapse_minus_diagnosis = relapse - diagnosis)
  p_wilcox <- if (nrow(wide) >= 2) stats::wilcox.test(wide$relapse, wide$diagnosis, paired = TRUE, exact = FALSE)$p.value else NA_real_
  p_t <- if (nrow(wide) >= 2) stats::t.test(wide$relapse, wide$diagnosis, paired = TRUE)$p.value else NA_real_
  p_sign <- if (nrow(wide) >= 1) stats::binom.test(sum(wide$delta_relapse_minus_diagnosis > 0, na.rm = TRUE), sum(!is.na(wide$delta_relapse_minus_diagnosis)), p = 0.5)$p.value else NA_real_
  list(
    wide = wide,
    stats = data.frame(
      metric = metric,
      n_pairs = nrow(wide),
      mean_diagnosis = mean(wide$diagnosis, na.rm = TRUE),
      mean_relapse = mean(wide$relapse, na.rm = TRUE),
      mean_delta_relapse_minus_diagnosis = mean(wide$delta_relapse_minus_diagnosis, na.rm = TRUE),
      pairs_increased = sum(wide$delta_relapse_minus_diagnosis > 0, na.rm = TRUE),
      paired_wilcox_p = p_wilcox,
      paired_t_p = p_t,
      sign_test_p = p_sign
    )
  )
}
metrics <- c("mean_STAT5B", "pct_STAT5B_detected", "STAT5B_pseudobulk_logCPM", "mean_proteasome_core", "mean_mitophagy_core")
comparison_list <- lapply(metrics, function(m) paired_compare(sample_summary, m))
paired_wide <- dplyr::bind_rows(lapply(seq_along(metrics), function(i) cbind(metric = metrics[i], comparison_list[[i]]$wide)))
paired_stats <- dplyr::bind_rows(lapply(comparison_list, `[[`, "stats"))
readr::write_csv(paired_wide, file.path(paths$results, "paired_diagnosis_relapse_metric_values.csv"))
readr::write_csv(paired_stats, file.path(paths$results, "paired_diagnosis_relapse_statistics.csv"))

plot_metric <- function(metric, ylab) {
  ggplot2::ggplot(sample_summary, ggplot2::aes(x = condition_simple, y = .data[[metric]], group = patient_pair, color = patient_pair)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_x_discrete(limits = c("diagnosis", "relapse")) +
    ggplot2::labs(x = NULL, y = ylab, color = "Pair") +
    theme_pub(11)
}
png(file.path(paths$figures, "paired_STAT5B_pathway_plots.png"), width = 1800, height = 1200, res = 150)
print(
  plot_metric("mean_STAT5B", "Mean STAT5B log-normalized expression") +
    plot_metric("STAT5B_pseudobulk_logCPM", "STAT5B pseudobulk logCPM") +
    plot_metric("mean_proteasome_core", "Mean proteasome score") +
    plot_metric("mean_mitophagy_core", "Mean mitophagy score")
)
grDevices::dev.off()

png(file.path(paths$figures, "STAT5B_feature_by_condition.png"), width = 1800, height = 750, res = 150)
print(Seurat::FeaturePlot(obj, features = "STAT5B", reduction = "umap.harmony", split.by = "condition_simple", raster = FALSE, ncol = 2))
grDevices::dev.off()
png(file.path(paths$figures, "STAT5B_violin_by_sample.png"), width = 1800, height = 750, res = 150)
print(Seurat::VlnPlot(obj, features = "STAT5B", group.by = "sample_id", pt.size = 0))
grDevices::dev.off()

saveRDS(obj, file.path(paths$objects, "GSE262271_seurat_qc_harmony_annotated.rds"))

report <- c(
  "# GSE262271 T-ALL Diagnosis-Relapse scRNA-seq Analysis Report",
  "",
  "## Dataset",
  "- GEO: GSE262271.",
  "- Source paper: Kypraios et al., Cancers 2024, DOI 10.3390/cancers16091667.",
  "- Three paired pediatric T-ALL diagnosis-relapse samples.",
  "- Main statistics are paired sample-level summaries across N=3 pairs, not pooled-cell p values.",
  "",
  "## Annotation Strategy",
  "- Manual marker-panel annotation only; no SingleR.",
  "- Annotation was checked using cluster-level DotPlot marker panels and Seurat FindAllMarkers top markers.",
  "- Broad downstream label is intentionally conservative: T_lineage_candidate_malignant versus B_cell_reference, Myeloid_reference, Erythroid_reference, pDC_reference, and Cytotoxic_T_or_NK_reference_like.",
  "- Detailed labels such as Cycling_TALL_blast_like, Quiescent_CD44_high_TALL_like, Immature_HSC_like_TALL, and ER_stress_TALL_like are retained as review aids, but downstream comparisons use the broad T-lineage candidate malignant compartment.",
  "",
  "## Outputs",
  "- QC: results/qc_filter_summary.csv",
  "- Annotation review: results/annotation_marker_panels_used.csv, results/cluster_top20_markers_res0.6_fast_review.csv, results/cluster_annotation_res0.6.csv, and figures/annotation_marker_dotplot_by_*.png",
  "- Candidate malignant-like cells: results/candidate_malignant_like_cells_metadata.csv",
  "- STAT5B/pathways: results/sample_level_candidate_malignant_like_summary.csv and results/paired_diagnosis_relapse_statistics.csv",
  "- inferCNV inputs and runner: infercnv_inputs/infercnv_raw_counts_matrix.rds, infercnv_inputs/infercnv_cell_annotations.txt, infercnv_inputs/infercnv_gene_order_from_orgHsEgDb.txt, and infercnv_inputs/run_infercnv_if_installed.R",
  "",
  "## Interpretation",
  "This dataset supports relapse-state association in T-ALL. It does not directly test bortezomib treatment response. Run infercnv_inputs/run_infercnv_if_installed.R and gse262271_infercnv_burden_summary.R after the downstream script to produce inferCNV heatmaps and CNV-burden summaries."
)
writeLines(report, file.path(paths$results, "ANALYSIS_REPORT.md"))
capture.output(sessionInfo(), file = file.path(paths$logs, paste0("sessionInfo_downstream_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")))
log_message("Downstream analysis complete")
