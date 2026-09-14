options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args) >= 1) args[[1]] else "outputs/GSE262271_TALL_relapse_scRNAseq_final", winslash = "/", mustWork = TRUE)
bulk_root <- normalizePath(if (length(args) >= 2) args[[2]] else Sys.getenv("GSE218858_BULK_ROOT", unset = "GSE218858_R_DESeq2_final_deliverables"), winslash = "/", mustWork = TRUE)
lib <- Sys.getenv("R_LIBS_USER", unset = "C:/Users/User/Documents/R/win-library/4.4")
if (dir.exists(lib)) .libPaths(unique(c(normalizePath(lib, winslash = "/"), .libPaths())))
for (p in c("dplyr", "ggplot2")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)

csv_dir <- file.path(root, "csv", "potential_malignant_T_cells")
png_dir <- file.path(root, "png", "potential_malignant_T_cells")
sc <- read.csv(file.path(csv_dir, "KEGG_GSEA_all_pathways_relapse_vs_diagnosis.csv"), check.names = FALSE)
bulk <- read.csv(file.path(bulk_root, "csv_results", "gsea", "R_fgsea_all_pathways_Hallmark_KEGG_GO_BP.csv"), check.names = FALSE)
bulk <- bulk[bulk$collection == "KEGG_2019_Mouse_Enrichr_mapped", , drop = FALSE]

map <- data.frame(
  display_pathway = c("Oxidative phosphorylation", "Mitophagy", "Proteasome", "Ribosome", "Ribosome biogenesis", "Spliceosome", "DNA replication"),
  sc_name = c("Oxidative phosphorylation", "Mitophagy - animal", "Proteasome", "Ribosome", "Ribosome biogenesis in eukaryotes", "Spliceosome", "DNA replication"),
  bulk_name = c("Oxidative phosphorylation", "Mitophagy", "Proteasome", "Ribosome", "Ribosome biogenesis in eukaryotes", "Spliceosome", "DNA replication"),
  stringsAsFactors = FALSE
)

split_genes <- function(x) {
  if (length(x) != 1 || is.na(x) || !nzchar(x)) return(character())
  unique(toupper(strsplit(x, ";", fixed = TRUE)[[1]]))
}

rows <- list()
for (i in seq_len(nrow(map))) {
  s <- sc[sc$pathway_name == map$sc_name[[i]], , drop = FALSE]
  b <- bulk[bulk$pathway == map$bulk_name[[i]], , drop = FALSE]
  sg <- split_genes(s$leadingEdge[[1]])
  bg <- split_genes(b$leadingEdge[[1]])
  overlap <- intersect(sg, bg)
  rows[[length(rows) + 1L]] <- data.frame(
    pathway = map$display_pathway[[i]],
    dataset = "Relapse scRNA-seq",
    contrast = "Relapse vs diagnosis; paired pseudobulk",
    NES = s$NES[[1]], FDR = s$padj[[1]], direction = s$direction[[1]],
    leading_edge_n = length(sg), cross_dataset_leading_edge_overlap_n = length(overlap),
    cross_dataset_leading_edge_overlap = paste(overlap, collapse = ";"),
    stringsAsFactors = FALSE
  )
  rows[[length(rows) + 1L]] <- data.frame(
    pathway = map$display_pathway[[i]],
    dataset = "STAT5B N642H bulk RNA-seq",
    contrast = "STAT5B N642H; Rag2-/- DN vs Rag2-/- DN",
    NES = b$NES[[1]], FDR = b$padj[[1]], direction = ifelse(b$NES[[1]] > 0, "STAT5B_N642H_up", "STAT5B_N642H_down"),
    leading_edge_n = length(bg), cross_dataset_leading_edge_overlap_n = length(overlap),
    cross_dataset_leading_edge_overlap = paste(overlap, collapse = ";"),
    stringsAsFactors = FALSE
  )
}
out <- do.call(rbind, rows)
out$significance <- ifelse(out$FDR < 0.05, "FDR < 0.05", ifelse(out$FDR < 0.10, "FDR < 0.10", "Not significant"))
write.csv(
  out,
  file.path(csv_dir, "cross_dataset_STAT5B_N642H_vs_relapse_GSEA_key_pathways.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# Gene-level audit for the four non-mitochondrial pathways emphasized by the user.
deg <- read.csv(file.path(csv_dir, "pseudobulk_DEG_relapse_vs_diagnosis_all_genes.csv"), check.names = FALSE)
stat5_sets <- read.csv(file.path(csv_dir, "STAT5_downstream_gene_sets_used.csv"), check.names = FALSE)
mito_le <- split_genes(sc$leadingEdge[sc$pathway_name == "Mitophagy - animal"][[1]])
prot_le <- split_genes(sc$leadingEdge[sc$pathway_name == "Proteasome"][[1]])
stat5_genes <- unique(toupper(stat5_sets$gene))
focus <- map[map$display_pathway %in% c("Ribosome", "Ribosome biogenesis", "Spliceosome", "DNA replication"), , drop = FALSE]
gene_rows <- list()
for (i in seq_len(nrow(focus))) {
  s <- sc[sc$pathway_name == focus$sc_name[[i]], , drop = FALSE]
  b <- bulk[bulk$pathway == focus$bulk_name[[i]], , drop = FALSE]
  sg <- split_genes(s$leadingEdge[[1]])
  bg <- split_genes(b$leadingEdge[[1]])
  gene_rows[[i]] <- data.frame(
    pathway = focus$display_pathway[[i]],
    gene = sg,
    shared_with_bulk_same_pathway = sg %in% bg,
    in_sc_mitophagy_leading_edge = sg %in% mito_le,
    in_sc_proteasome_leading_edge = sg %in% prot_le,
    in_curated_STAT5_target_sets = sg %in% stat5_genes,
    stringsAsFactors = FALSE
  )
}
gene_audit <- do.call(rbind, gene_rows)
deg2 <- deg[, c("gene", "logFC", "PValue", "FDR", "rank_stat")]
deg2$gene <- toupper(deg2$gene)
gene_audit <- dplyr::left_join(gene_audit, deg2, by = "gene") |>
  dplyr::arrange(pathway, dplyr::desc(rank_stat))
write.csv(
  gene_audit,
  file.path(csv_dir, "ribosome_spliceosome_DNAreplication_leading_edge_gene_audit.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

out$pathway <- factor(out$pathway, levels = rev(map$display_pathway))
out$dataset <- factor(out$dataset, levels = c("Relapse scRNA-seq", "STAT5B N642H bulk RNA-seq"))
p <- ggplot2::ggplot(out, ggplot2::aes(x = NES, y = pathway, color = significance, shape = dataset)) +
  ggplot2::geom_vline(xintercept = 0, color = "grey70", linewidth = 0.45) +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = NES, yend = pathway), color = "grey80", linewidth = 0.45, position = ggplot2::position_dodge(width = 0.55)) +
  ggplot2::geom_point(size = 3.0, stroke = 0.8, position = ggplot2::position_dodge(width = 0.55)) +
  ggplot2::scale_color_manual(values = c("FDR < 0.05" = "#B2433F", "FDR < 0.10" = "#D89166", "Not significant" = "grey65")) +
  ggplot2::scale_shape_manual(values = c("Relapse scRNA-seq" = 16, "STAT5B N642H bulk RNA-seq" = 17)) +
  ggplot2::labs(
    title = "Pathway comparison across relapse T-ALL and STAT5B N642H models",
    subtitle = "Positive NES indicates enrichment in relapse or STAT5B N642H cells",
    x = "Normalized enrichment score (NES)", y = NULL, color = NULL, shape = NULL,
    caption = "GSEA results are shown within each dataset; FDR values are not pooled across datasets."
  ) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    axis.text = ggplot2::element_text(color = "black"),
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "bottom"
  )
ggplot2::ggsave(
  file.path(png_dir, "cross_dataset_STAT5B_N642H_vs_relapse_GSEA_key_pathways.png"),
  p, width = 8.2, height = 5.4, dpi = 400, bg = "white"
)

print(out[, c("pathway", "dataset", "NES", "FDR", "significance", "cross_dataset_leading_edge_overlap_n")])
