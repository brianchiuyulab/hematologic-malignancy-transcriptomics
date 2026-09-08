options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args) >= 1) args[[1]] else "outputs/GSE262271_TALL_relapse_scRNAseq_final", winslash = "/", mustWork = TRUE)
rds <- normalizePath(if (length(args) >= 2) args[[2]] else "04_R_Objects/GSE262271_seurat_qc_harmony_annotated.rds", winslash = "/", mustWork = TRUE)
lib <- Sys.getenv("R_LIBS_USER", unset = "C:/Users/User/Documents/R/win-library/4.4")
if (dir.exists(lib)) .libPaths(unique(c(normalizePath(lib, winslash = "/"), .libPaths())))
for (p in c("Seurat", "SeuratObject", "Matrix", "edgeR", "dplyr")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)

csv_dir <- file.path(root, "csv", "potential_malignant_T_cells")
gsea <- read.csv(file.path(csv_dir, "KEGG_GSEA_all_pathways_relapse_vs_diagnosis.csv"), check.names = FALSE)
sel <- read.csv(file.path(csv_dir, "potential_malignant_T_cells_metadata.csv"), check.names = FALSE)
paths <- c("Mitophagy - animal", "Autophagy - other", "Autophagy - animal")
gsea <- gsea[gsea$pathway_name %in% paths, , drop = FALSE]

obj <- readRDS(rds)
obj <- subset(obj, cells = intersect(sel$cell, colnames(obj)))
obj <- tryCatch(SeuratObject::JoinLayers(obj, assay = "RNA"), error = function(e) obj)
counts <- tryCatch(SeuratObject::LayerData(obj[["RNA"]], layer = "counts"), error = function(e) Seurat::GetAssayData(obj, assay = "RNA", slot = "counts"))
meta <- obj@meta.data
meta$cell <- rownames(meta)
meta <- meta[match(colnames(counts), meta$cell), , drop = FALSE]
samples <- unique(meta$sample_id[order(meta$patient_pair, meta$condition_simple)])
pb <- do.call(cbind, lapply(samples, function(s) Matrix::rowSums(counts[, meta$sample_id == s, drop = FALSE])))
rownames(pb) <- rownames(counts); colnames(pb) <- samples
sm <- do.call(rbind, lapply(samples, function(s) {
  z <- meta[meta$sample_id == s, , drop = FALSE][1, ]
  data.frame(sample_id = s, patient_pair = z$patient_pair, condition_simple = z$condition_simple)
}))
y <- edgeR::calcNormFactors(edgeR::DGEList(pb))
lcpm <- edgeR::cpm(y, log = TRUE, prior.count = 2)

rows <- list(); k <- 0L
for (i in seq_len(nrow(gsea))) {
  genes <- intersect(strsplit(gsea$leadingEdge[[i]], ";", fixed = TRUE)[[1]], rownames(lcpm))
  for (gene in genes) {
    delta <- setNames(numeric(3), c("P1", "P2", "P3"))
    for (patient in names(delta)) {
      d <- sm$sample_id[sm$patient_pair == patient & sm$condition_simple == "diagnosis"]
      r <- sm$sample_id[sm$patient_pair == patient & sm$condition_simple == "relapse"]
      delta[[patient]] <- lcpm[gene, r] - lcpm[gene, d]
    }
    k <- k + 1L
    rows[[k]] <- data.frame(
      pathway_name = gsea$pathway_name[[i]], gene = gene,
      P1_delta = delta[["P1"]], P2_delta = delta[["P2"]], P3_delta = delta[["P3"]],
      n_pairs_up = sum(delta > 0), n_pairs_down = sum(delta < 0),
      direction_consistency = if (sum(delta > 0) == 3) "up in 3/3" else if (sum(delta > 0) == 2) "up in 2/3" else if (sum(delta < 0) == 3) "down in 3/3" else if (sum(delta < 0) == 2) "down in 2/3" else "mixed",
      stringsAsFactors = FALSE
    )
  }
}
out <- do.call(rbind, rows)
write.csv(out, file.path(csv_dir, "mitophagy_autophagy_leading_edge_paired_consistency.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(
  out[out$pathway_name == "Mitophagy - animal", , drop = FALSE],
  file.path(csv_dir, "mitophagy_leading_edge_paired_directions.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
summary <- out |>
  dplyr::group_by(pathway_name) |>
  dplyr::summarise(
    leading_edge_n = dplyr::n(), up_3_of_3 = sum(n_pairs_up == 3), up_2_of_3 = sum(n_pairs_up == 2),
    down_2_of_3 = sum(n_pairs_down == 2), down_3_of_3 = sum(n_pairs_down == 3), .groups = "drop"
  ) |>
  dplyr::left_join(gsea[, c("pathway_name", "NES", "pval", "padj", "direction")], by = "pathway_name")
write.csv(summary, file.path(csv_dir, "mitophagy_autophagy_paired_consistency_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")
print(summary)
