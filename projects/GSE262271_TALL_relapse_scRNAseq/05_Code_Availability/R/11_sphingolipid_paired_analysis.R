options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(args) >= 1) args[[1]] else "outputs/GSE262271_TALL_relapse_scRNAseq_final"
seurat_file <- if (length(args) >= 2) args[[2]] else file.path(
  "C:/Users/User/Desktop/Public dataset analysis/GEO dataset/2026/single cell/BTZ 相關",
  "GSE262271_TALL_relapse_scRNAseq/processed_objects/GSE262271_seurat_qc_harmony_annotated.rds"
)
output_root <- normalizePath(output_root, winslash = "/", mustWork = TRUE)
seurat_file <- normalizePath(seurat_file, winslash = "/", mustWork = TRUE)

user_lib <- Sys.getenv("R_LIBS_USER", unset = "C:/Users/User/Documents/R/win-library/4.4")
if (dir.exists(user_lib)) .libPaths(unique(c(normalizePath(user_lib, winslash = "/"), .libPaths())))
required <- c("Seurat", "SeuratObject", "Matrix", "edgeR", "dplyr", "ggplot2", "patchwork")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "), call. = FALSE)

csv_dir <- file.path(output_root, "csv", "potential_malignant_T_cells")
png_dir <- file.path(output_root, "png", "potential_malignant_T_cells")
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)

gsea <- utils::read.csv(file.path(csv_dir, "KEGG_GSEA_all_pathways_relapse_vs_diagnosis.csv"), check.names = FALSE)
deg <- utils::read.csv(file.path(csv_dir, "pseudobulk_DEG_relapse_vs_diagnosis_all_genes.csv"), check.names = FALSE)
selected_cells <- utils::read.csv(file.path(csv_dir, "potential_malignant_T_cells_metadata.csv"), check.names = FALSE)

target_pathways <- c(
  "Sphingolipid signaling pathway",
  "Sphingolipid metabolism",
  "Glycosphingolipid biosynthesis - lacto and neolacto series",
  "Glycosphingolipid biosynthesis - globo and isoglobo series",
  "Glycosphingolipid biosynthesis - ganglio series"
)
pathway_summary <- gsea[gsea$pathway_name %in% target_pathways, c(
  "pathway_id", "pathway_name", "NES", "pval", "padj", "size", "direction", "leadingEdge", "linked_genes"
), drop = FALSE]
pathway_summary <- pathway_summary[match(target_pathways, pathway_summary$pathway_name, nomatch = 0), , drop = FALSE]
names(pathway_summary)[names(pathway_summary) == "pval"] <- "GSEA_PValue"
names(pathway_summary)[names(pathway_summary) == "padj"] <- "GSEA_FDR"
utils::write.csv(
  pathway_summary,
  file.path(csv_dir, "sphingolipid_KEGG_GSEA_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

message("Reading annotated Seurat object...")
obj <- readRDS(seurat_file)
keep_cells <- intersect(selected_cells$cell, colnames(obj))
if (!length(keep_cells)) stop("No selected potential malignant T cells were found in the Seurat object.")
obj <- subset(obj, cells = keep_cells)
obj <- tryCatch(SeuratObject::JoinLayers(obj, assay = "RNA"), error = function(e) obj)
counts <- tryCatch(
  SeuratObject::LayerData(obj[["RNA"]], layer = "counts"),
  error = function(e) Seurat::GetAssayData(obj, assay = "RNA", slot = "counts")
)

cell_meta <- obj@meta.data
cell_meta$cell <- rownames(cell_meta)
cell_meta <- cell_meta[match(colnames(counts), cell_meta$cell), , drop = FALSE]
sample_levels <- c("M104", "M127", "M105", "M128", "M106", "M129")
sample_levels <- sample_levels[sample_levels %in% unique(cell_meta$sample_id)]
if (length(sample_levels) != 6) {
  sample_levels <- unique(cell_meta$sample_id[order(cell_meta$patient_pair, cell_meta$condition_simple)])
}

message("Aggregating counts to six patient-level pseudobulk samples...")
pb_counts <- do.call(cbind, lapply(sample_levels, function(s) {
  Matrix::rowSums(counts[, cell_meta$sample_id == s, drop = FALSE])
}))
rownames(pb_counts) <- rownames(counts)
colnames(pb_counts) <- sample_levels

sample_meta <- do.call(rbind, lapply(sample_levels, function(s) {
  z <- cell_meta[cell_meta$sample_id == s, , drop = FALSE][1, , drop = FALSE]
  data.frame(
    sample_id = s,
    patient_pair = as.character(z$patient_pair),
    condition_simple = as.character(z$condition_simple),
    n_cells = sum(cell_meta$sample_id == s),
    stringsAsFactors = FALSE
  )
}))

y <- edgeR::DGEList(counts = pb_counts)
y <- edgeR::calcNormFactors(y)
logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 2)

split_genes <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character())
  unique(strsplit(x, ";", fixed = TRUE)[[1]])
}
all_membership <- do.call(rbind, lapply(seq_len(nrow(pathway_summary)), function(i) {
  data.frame(
    pathway_name = pathway_summary$pathway_name[[i]],
    gene = split_genes(pathway_summary$linked_genes[[i]]),
    is_leading_edge = FALSE,
    stringsAsFactors = FALSE
  )
}))
leading_membership <- do.call(rbind, lapply(seq_len(nrow(pathway_summary)), function(i) {
  data.frame(
    pathway_name = pathway_summary$pathway_name[[i]],
    gene = split_genes(pathway_summary$leadingEdge[[i]]),
    is_leading_edge = TRUE,
    stringsAsFactors = FALSE
  )
}))
membership <- dplyr::bind_rows(all_membership, leading_membership) |>
  dplyr::group_by(pathway_name, gene) |>
  dplyr::summarise(is_leading_edge = any(is_leading_edge), .groups = "drop") |>
  dplyr::filter(gene %in% rownames(logcpm))

target_genes <- sort(unique(membership$gene))
expr_long <- do.call(rbind, lapply(target_genes, function(g) {
  data.frame(
    gene = g,
    sample_meta,
    pseudobulk_logCPM = as.numeric(logcpm[g, sample_meta$sample_id]),
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(
  expr_long,
  file.path(csv_dir, "sphingolipid_gene_sample_pseudobulk_logCPM.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

paired_delta <- do.call(rbind, lapply(split(expr_long, list(expr_long$gene, expr_long$patient_pair), drop = TRUE), function(z) {
  d <- z$pseudobulk_logCPM[z$condition_simple == "diagnosis"]
  r <- z$pseudobulk_logCPM[z$condition_simple == "relapse"]
  if (length(d) != 1 || length(r) != 1) return(NULL)
  data.frame(
    gene = z$gene[[1]],
    patient_pair = z$patient_pair[[1]],
    relapse_minus_diagnosis_logCPM = r - d,
    direction = ifelse(r > d, "up", ifelse(r < d, "down", "unchanged")),
    stringsAsFactors = FALSE
  )
}))

gene_summary <- paired_delta |>
  dplyr::group_by(gene) |>
  dplyr::summarise(
    P1_delta = relapse_minus_diagnosis_logCPM[match("P1", patient_pair)],
    P2_delta = relapse_minus_diagnosis_logCPM[match("P2", patient_pair)],
    P3_delta = relapse_minus_diagnosis_logCPM[match("P3", patient_pair)],
    n_pairs_up = sum(relapse_minus_diagnosis_logCPM > 0),
    n_pairs_down = sum(relapse_minus_diagnosis_logCPM < 0),
    mean_pair_delta = mean(relapse_minus_diagnosis_logCPM),
    median_pair_delta = stats::median(relapse_minus_diagnosis_logCPM),
    .groups = "drop"
  ) |>
  dplyr::left_join(deg |> dplyr::select(gene, logFC, PValue, FDR, rank_stat), by = "gene") |>
  dplyr::left_join(
    membership |>
      dplyr::group_by(gene) |>
      dplyr::summarise(
        pathway_membership = paste(pathway_name, collapse = "; "),
        leading_edge_pathways = paste(pathway_name[is_leading_edge], collapse = "; "),
        any_leading_edge = any(is_leading_edge),
        .groups = "drop"
      ),
    by = "gene"
  ) |>
  dplyr::mutate(
    paired_direction_consistency = dplyr::case_when(
      n_pairs_up == 3 ~ "up in 3/3",
      n_pairs_down == 3 ~ "down in 3/3",
      n_pairs_up == 2 ~ "up in 2/3",
      n_pairs_down == 2 ~ "down in 2/3",
      TRUE ~ "mixed"
    ),
    gene_level_FDR_lt_0.05 = FDR < 0.05
  ) |>
  dplyr::arrange(dplyr::desc(any_leading_edge), dplyr::desc(logFC))

utils::write.csv(
  gene_summary,
  file.path(csv_dir, "sphingolipid_gene_paired_directions.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

leading_summary <- gene_summary |>
  dplyr::filter(any_leading_edge) |>
  dplyr::arrange(dplyr::desc(logFC))
utils::write.csv(
  leading_summary,
  file.path(csv_dir, "sphingolipid_leading_edge_gene_paired_directions.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# Publication-style paired-direction plot: each tile is one patient's R-D change.
plot_genes <- leading_summary$gene
plot_delta <- paired_delta[paired_delta$gene %in% plot_genes, , drop = FALSE]
plot_delta$gene <- factor(plot_delta$gene, levels = rev(plot_genes))
plot_delta$patient_pair <- factor(plot_delta$patient_pair, levels = c("P1", "P2", "P3"))
lim <- stats::quantile(abs(plot_delta$relapse_minus_diagnosis_logCPM), 0.95, na.rm = TRUE)
lim <- max(1, min(4, as.numeric(lim)))
plot_delta$delta_plot <- pmax(-lim, pmin(lim, plot_delta$relapse_minus_diagnosis_logCPM))

p_heat <- ggplot2::ggplot(plot_delta, ggplot2::aes(patient_pair, gene, fill = delta_plot)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.25) +
  ggplot2::scale_fill_gradient2(
    low = "#3C6E9E", mid = "white", high = "#B2433F", midpoint = 0,
    limits = c(-lim, lim), name = "R-D\nlogCPM"
  ) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(
    axis.text = ggplot2::element_text(color = "black"),
    axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    legend.position = "bottom",
    plot.margin = ggplot2::margin(4, 2, 4, 4)
  )

side <- leading_summary
side$gene <- factor(side$gene, levels = rev(plot_genes))
side$consistency <- factor(
  side$paired_direction_consistency,
  levels = c("up in 3/3", "up in 2/3", "down in 2/3", "down in 3/3", "mixed")
)
p_fc <- ggplot2::ggplot(side, ggplot2::aes(logFC, gene, color = consistency)) +
  ggplot2::geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = logFC, yend = gene), linewidth = 0.45) +
  ggplot2::geom_point(size = 1.8) +
  ggplot2::scale_color_manual(values = c(
    "up in 3/3" = "#9D2727", "up in 2/3" = "#D7785D",
    "down in 2/3" = "#6D93B3", "down in 3/3" = "#285D88", "mixed" = "grey55"
  ), drop = FALSE) +
  ggplot2::labs(x = "Paired model log2FC", y = NULL, color = "Direction") +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.line.y = ggplot2::element_blank(),
    legend.position = "bottom",
    plot.margin = ggplot2::margin(4, 4, 4, 2)
  )

combined <- (p_heat | p_fc) +
  patchwork::plot_layout(widths = c(1, 1.5), guides = "collect") +
  patchwork::plot_annotation(
    title = "Sphingolipid pathway leading-edge genes",
    subtitle = "Potential malignant T cells; paired diagnosis-relapse samples (N = 3)",
    caption = "Tiles show within-patient relapse minus diagnosis changes. No individual gene passed FDR < 0.05.",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 9),
      plot.caption = ggplot2::element_text(size = 8, hjust = 0)
    )
  ) & ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(
  file.path(png_dir, "sphingolipid_leading_edge_paired_directions.png"),
  combined,
  width = 7.2,
  height = max(6.5, 0.19 * length(plot_genes) + 2.4),
  dpi = 400,
  bg = "white"
)

cat("Potential malignant T cells:", length(keep_cells), "\n")
cat("Sphingolipid leading-edge genes:", nrow(leading_summary), "\n")
cat("Leading-edge genes up in 3/3 pairs:", sum(leading_summary$n_pairs_up == 3), "\n")
cat("Leading-edge genes down in 3/3 pairs:", sum(leading_summary$n_pairs_down == 3), "\n")
cat("Gene-level FDR < 0.05:", sum(leading_summary$FDR < 0.05, na.rm = TRUE), "\n")
