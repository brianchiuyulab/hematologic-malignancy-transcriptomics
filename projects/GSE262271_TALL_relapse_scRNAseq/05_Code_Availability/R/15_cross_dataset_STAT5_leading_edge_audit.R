options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args) >= 1) args[[1]] else "outputs/GSE262271_TALL_relapse_scRNAseq_final", winslash = "/", mustWork = TRUE)
bulk_root <- normalizePath(if (length(args) >= 2) args[[2]] else Sys.getenv("GSE218858_BULK_ROOT", unset = "GSE218858_R_DESeq2_final_deliverables"), winslash = "/", mustWork = TRUE)
lib <- Sys.getenv("R_LIBS_USER", unset = "C:/Users/User/Documents/R/win-library/4.4")
if (dir.exists(lib)) .libPaths(unique(c(normalizePath(lib, winslash = "/"), .libPaths())))
if (!requireNamespace("dplyr", quietly = TRUE)) stop("Missing package: dplyr")

csv_dir <- file.path(root, "csv", "potential_malignant_T_cells")
sc <- read.csv(file.path(csv_dir, "KEGG_GSEA_all_pathways_relapse_vs_diagnosis.csv"), check.names = FALSE)
bulk <- read.csv(file.path(bulk_root, "csv_results", "gsea", "R_fgsea_all_pathways_Hallmark_KEGG_GO_BP.csv"), check.names = FALSE)
bulk <- bulk[bulk$collection == "KEGG_2019_Mouse_Enrichr_mapped", , drop = FALSE]
shared <- read.csv(file.path(csv_dir, "cross_dataset_KEGG_FDR0.10_same_direction_intersection.csv"), check.names = FALSE)
stat5_sets <- read.csv(file.path(csv_dir, "STAT5_downstream_gene_sets_used.csv"), check.names = FALSE)

canon <- function(x) {
  x <- trimws(x)
  x[x == "Mitophagy - animal"] <- "Mitophagy"
  x[x == "Autophagy - animal"] <- "Autophagy"
  x
}
split_genes <- function(x) {
  if (length(x) != 1 || is.na(x) || !nzchar(x)) return(character())
  unique(toupper(strsplit(x, ";", fixed = TRUE)[[1]]))
}

sc$canonical <- canon(sc$pathway_name)
bulk$canonical <- canon(bulk$pathway)
curated_stat5 <- unique(toupper(stat5_sets$gene))
stat5_axis <- unique(c(
  "STAT5A", "STAT5B", "JAK1", "JAK2", "JAK3", "TYK2",
  "IL2RA", "IL2RB", "IL2RG", "IL4R", "IL7R", "IL9R", "IL15RA", "IL21R",
  "IFNGR1", "IFNGR2", "OSM", "LIF", "CSF2", "CSF2RB",
  "CISH", "SOCS1", "SOCS2", "SOCS3", "SOCS6", "PIM1", "PIM2", "PIM3",
  "BCL2", "BCL2L1", "MCL1", "MYC", "CCND1", "CCND2", "CCND3"
))

focus_extra <- c("VEGF signaling pathway")
pathways <- unique(c(shared$pathway, focus_extra))
summary_rows <- list()
gene_rows <- list()
for (p in pathways) {
  s <- sc[sc$canonical == p, , drop = FALSE]
  b <- bulk[bulk$canonical == p, , drop = FALSE]
  if (!nrow(s) || !nrow(b)) next
  s <- s[which.min(s$padj), , drop = FALSE]
  b <- b[which.min(b$padj), , drop = FALSE]
  sg <- split_genes(s$leadingEdge[[1]])
  bg <- split_genes(b$leadingEdge[[1]])
  ov <- intersect(sg, bg)
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    pathway = p,
    in_same_direction_FDR0.10_intersection = p %in% shared$pathway,
    sc_NES = s$NES[[1]], sc_FDR = s$padj[[1]],
    bulk_NES = b$NES[[1]], bulk_FDR = b$padj[[1]],
    sc_leading_edge_n = length(sg), bulk_leading_edge_n = length(bg), shared_leading_edge_n = length(ov),
    shared_leading_edge_genes = paste(ov, collapse = ";"),
    shared_STAT5_axis_genes = paste(intersect(ov, stat5_axis), collapse = ";"),
    sc_STAT5_axis_genes = paste(intersect(sg, stat5_axis), collapse = ";"),
    bulk_STAT5_axis_genes = paste(intersect(bg, stat5_axis), collapse = ";"),
    shared_curated_STAT5_targets = paste(intersect(ov, curated_stat5), collapse = ";"),
    stringsAsFactors = FALSE
  )
  union_genes <- union(sg, bg)
  gene_rows[[length(gene_rows) + 1L]] <- data.frame(
    pathway = p, gene = union_genes,
    in_sc_leading_edge = union_genes %in% sg,
    in_bulk_leading_edge = union_genes %in% bg,
    shared_leading_edge = union_genes %in% ov,
    STAT5_axis_gene = union_genes %in% stat5_axis,
    curated_STAT5_target = union_genes %in% curated_stat5,
    stringsAsFactors = FALSE
  )
}

summary_out <- do.call(rbind, summary_rows) |>
  dplyr::arrange(dplyr::desc(in_same_direction_FDR0.10_intersection), pmax(sc_FDR, bulk_FDR), pathway)
gene_out <- do.call(rbind, gene_rows) |>
  dplyr::arrange(pathway, dplyr::desc(shared_leading_edge), gene)
write.csv(summary_out, file.path(csv_dir, "cross_dataset_STAT5_pathway_leading_edge_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(gene_out, file.path(csv_dir, "cross_dataset_STAT5_pathway_leading_edge_gene_audit.csv"), row.names = FALSE, fileEncoding = "UTF-8")

recurrent <- gene_out |>
  dplyr::filter(shared_leading_edge) |>
  dplyr::group_by(gene) |>
  dplyr::summarise(
    shared_pathways_n = dplyr::n_distinct(pathway),
    pathways = paste(sort(unique(pathway)), collapse = ";"),
    STAT5_axis_gene = any(STAT5_axis_gene),
    curated_STAT5_target = any(curated_STAT5_target),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(shared_pathways_n), gene)
write.csv(recurrent, file.path(csv_dir, "cross_dataset_recurrent_shared_leading_edge_genes.csv"), row.names = FALSE, fileEncoding = "UTF-8")

print(summary_out[summary_out$pathway %in% c(
  "VEGF signaling pathway", "JAK-STAT signaling pathway", "Th1 and Th2 cell differentiation",
  "Th17 cell differentiation", "B cell receptor signaling pathway", "Rap1 signaling pathway",
  "T cell receptor signaling pathway", "HIF-1 signaling pathway", "Mitophagy"
), c("pathway", "sc_FDR", "bulk_FDR", "shared_leading_edge_n", "shared_STAT5_axis_genes")], row.names = FALSE)
cat("\nMost recurrent shared leading-edge genes:\n")
print(head(recurrent, 30), row.names = FALSE)
