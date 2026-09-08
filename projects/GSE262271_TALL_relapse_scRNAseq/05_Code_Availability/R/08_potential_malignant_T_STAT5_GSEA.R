options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = "") {
  hit <- match(flag, args)
  if (!is.na(hit) && length(args) >= hit + 1) args[[hit + 1]] else default
}

output_root <- normalizePath(
  get_arg("--output-root", "outputs/GSE262271_TALL_relapse_scRNAseq_final"),
  winslash = "/",
  mustWork = TRUE
)

user_lib <- Sys.getenv("R_LIBS_USER", unset = "C:/Users/User/Documents/R/win-library/4.4")
if (dir.exists(user_lib)) .libPaths(unique(c(normalizePath(user_lib, winslash = "/"), .libPaths())))

required <- c("dplyr", "ggplot2", "fgsea")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "), call. = FALSE)

csv_root <- file.path(output_root, "csv")
png_root <- file.path(output_root, "png")
final_csv <- file.path(csv_root, "potential_malignant_T_cells")
final_png <- file.path(png_root, "potential_malignant_T_cells")
dir.create(final_csv, recursive = TRUE, showWarnings = FALSE)
dir.create(final_png, recursive = TRUE, showWarnings = FALSE)

deg_file <- file.path(
  csv_root,
  "sensitivity_mitophagy_gsea",
  "DEG_pseudobulk_CNV_high_top25_candidate_relapse_vs_diagnosis.csv"
)
gene_set_file <- file.path(csv_root, "STAT5_target_gene_sets_used.csv")
count_file <- file.path(csv_root, "sensitivity_mitophagy_gsea", "sensitivity_compartment_cell_counts.csv")
kegg_file <- file.path(csv_root, "sensitivity_mitophagy_gsea", "sensitivity_kegg_all_pathways_by_compartment.csv")

stopifnot(file.exists(deg_file), file.exists(gene_set_file), file.exists(count_file), file.exists(kegg_file))

deg <- utils::read.csv(deg_file, check.names = FALSE)
deg <- deg[is.finite(deg$rank_stat) & !is.na(deg$gene) & nzchar(deg$gene), , drop = FALSE]
deg <- deg[order(-abs(deg$rank_stat)), , drop = FALSE]
deg <- deg[!duplicated(deg$gene), , drop = FALSE]
rank_stats <- deg$rank_stat
names(rank_stats) <- deg$gene
rank_stats <- sort(rank_stats, decreasing = TRUE)

gene_set_table <- utils::read.csv(gene_set_file, check.names = FALSE)
gene_sets <- split(gene_set_table$gene, gene_set_table$target_set)
gene_sets <- lapply(gene_sets, function(x) intersect(unique(x), names(rank_stats)))
gene_sets <- gene_sets[lengths(gene_sets) >= 3]

set.seed(20260716)
gsea <- fgsea::fgseaMultilevel(
  pathways = gene_sets,
  stats = rank_stats,
  minSize = 3,
  maxSize = 500,
  eps = 0
)
gsea$leadingEdge <- vapply(gsea$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
gsea <- as.data.frame(gsea)
gsea$direction <- ifelse(gsea$NES > 0, "relapse_up", ifelse(gsea$NES < 0, "relapse_down", "flat"))
gsea$significant_FDR_0.05 <- gsea$padj < 0.05
gsea$analysis_population <- "potential malignant T cells (CNV-high top 25%)"
gsea$comparison <- "relapse_vs_diagnosis_paired_pseudobulk_N3"
gsea <- gsea[, c(
  "analysis_population", "comparison", "pathway", "size", "ES", "NES", "pval", "padj",
  "log2err", "direction", "significant_FDR_0.05", "leadingEdge"
)]
gsea <- gsea[order(gsea$padj, -abs(gsea$NES)), , drop = FALSE]
utils::write.csv(
  gsea,
  file.path(final_csv, "STAT5_downstream_gene_set_GSEA_relapse_vs_diagnosis.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

gene_set_table$analysis_population <- "potential malignant T cells (CNV-high top 25%)"
utils::write.csv(
  gene_set_table[, c("analysis_population", "target_set", "gene")],
  file.path(final_csv, "STAT5_downstream_gene_sets_used.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

deg$analysis_population <- "potential malignant T cells (CNV-high top 25%)"
utils::write.csv(
  deg[, c("analysis_population", "gene", "logFC", "logCPM", "F", "PValue", "FDR", "rank_stat")],
  file.path(final_csv, "pseudobulk_DEG_relapse_vs_diagnosis_all_genes.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

stat5_ab <- deg[deg$gene %in% c("STAT5A", "STAT5B"), c("analysis_population", "gene", "logFC", "logCPM", "F", "PValue", "FDR", "rank_stat")]
utils::write.csv(
  stat5_ab,
  file.path(final_csv, "STAT5A_STAT5B_pseudobulk_DEG.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

kegg <- utils::read.csv(kegg_file, check.names = FALSE)
kegg_top25 <- kegg[kegg$compartment == "CNV_high_top25_candidate", , drop = FALSE]
kegg_top25$analysis_population <- "potential malignant T cells (CNV-high top 25%)"
kegg_top25$compartment <- "potential_malignant_T_cells_CNV_high_top25"
kegg_top25 <- kegg_top25[, c("analysis_population", setdiff(names(kegg_top25), "analysis_population"))]
utils::write.csv(
  kegg_top25,
  file.path(final_csv, "KEGG_GSEA_all_pathways_relapse_vs_diagnosis.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

filter_top25_result <- function(input_name, output_name) {
  x <- utils::read.csv(file.path(csv_root, "sensitivity_mitophagy_gsea", input_name), check.names = FALSE)
  x <- x[x$compartment == "CNV_high_top25_candidate", , drop = FALSE]
  x$analysis_population <- "potential malignant T cells (CNV-high top 25%)"
  x$compartment <- "potential_malignant_T_cells_CNV_high_top25"
  x <- x[, c("analysis_population", setdiff(names(x), "analysis_population"))]
  utils::write.csv(x, file.path(final_csv, output_name), row.names = FALSE, fileEncoding = "UTF-8")
  x
}

kegg_focus_top25 <- filter_top25_result(
  "sensitivity_kegg_focus_pathways_by_compartment.csv",
  "KEGG_GSEA_focus_pathways_relapse_vs_diagnosis.csv"
)
go_focus_top25 <- filter_top25_result(
  "sensitivity_go_bp_focus_terms_by_compartment.csv",
  "GO_BP_GSEA_mito_autophagy_proteasome_focus.csv"
)
custom_mito_top25 <- filter_top25_result(
  "sensitivity_custom_mito_gene_sets_fgsea_by_compartment.csv",
  "custom_mitochondrial_gene_set_GSEA.csv"
)
custom_module_stats_top25 <- filter_top25_result(
  "sensitivity_custom_mito_module_paired_statistics_by_compartment.csv",
  "custom_mitochondrial_module_paired_statistics.csv"
)
custom_module_samples_top25 <- filter_top25_result(
  "sensitivity_custom_mito_module_sample_scores_by_compartment.csv",
  "custom_mitochondrial_module_sample_scores.csv"
)
mito_summary_top25 <- filter_top25_result(
  "sensitivity_summary_mitophagy_mito_proteasome.csv",
  "mitophagy_mitochondrial_proteasome_summary.csv"
)

heatmap_csv_root <- file.path(csv_root, "proteasome_mitophagy_heatmaps")
heatmap_sample <- utils::read.csv(
  file.path(heatmap_csv_root, "KEGG_proteasome_mitophagy_sample_pseudobulk_logCPM.csv"),
  check.names = FALSE
)
heatmap_sample <- heatmap_sample[heatmap_sample$compartment == "CNV_high_top25_candidate", , drop = FALSE]
heatmap_sample$analysis_population <- "potential malignant T cells (CNV-high top 25%)"
heatmap_sample$compartment <- "potential_malignant_T_cells_CNV_high_top25"
heatmap_sample <- heatmap_sample[, c("analysis_population", setdiff(names(heatmap_sample), "analysis_population"))]
utils::write.csv(
  heatmap_sample,
  file.path(final_csv, "proteasome_mitophagy_sample_pseudobulk_logCPM.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
heatmap_deg <- utils::read.csv(
  file.path(heatmap_csv_root, "KEGG_proteasome_mitophagy_DEG_by_compartment.csv"),
  check.names = FALSE
)
heatmap_deg <- heatmap_deg[heatmap_deg$compartment == "CNV_high_top25_candidate", , drop = FALSE]
heatmap_deg$analysis_population <- "potential malignant T cells (CNV-high top 25%)"
heatmap_deg$compartment <- "potential_malignant_T_cells_CNV_high_top25"
heatmap_deg <- heatmap_deg[, c("analysis_population", setdiff(names(heatmap_deg), "analysis_population"))]
utils::write.csv(
  heatmap_deg,
  file.path(final_csv, "proteasome_mitophagy_gene_DEG.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
source_heatmap_png <- file.path(
  png_root,
  "proteasome_mitophagy_heatmaps",
  "KEGG_proteasome_mitophagy_CNV_high_top25_sample_heatmap.png"
)
if (file.exists(source_heatmap_png)) {
  file.copy(
    source_heatmap_png,
    file.path(final_png, "proteasome_mitophagy_sample_heatmap.png"),
    overwrite = TRUE
  )
}

cnv_file <- file.path(csv_root, "infercnv_nonTref_cell_cnv_burden.csv")
cnv <- utils::read.csv(cnv_file, check.names = FALSE)
candidate_cnv <- cnv[
  cnv$broad_annotation == "T_lineage_candidate_malignant" & is.finite(cnv$infercnv_mean_abs_log2),
  ,
  drop = FALSE
]
cnv_q75 <- as.numeric(stats::quantile(candidate_cnv$infercnv_mean_abs_log2, 0.75, na.rm = TRUE))
potential_cells <- candidate_cnv[candidate_cnv$infercnv_mean_abs_log2 >= cnv_q75, , drop = FALSE]
potential_cells$analysis_population <- "potential malignant T cells (CNV-high top 25%)"
potential_cells$selection_threshold_q75 <- cnv_q75
potential_cells <- potential_cells[, c("analysis_population", "selection_threshold_q75", setdiff(names(potential_cells), c("analysis_population", "selection_threshold_q75")))]
utils::write.csv(
  potential_cells,
  file.path(final_csv, "potential_malignant_T_cells_metadata.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
selection_summary <- data.frame(
  analysis_population = "potential malignant T cells",
  operational_definition = "T-lineage candidate malignant-like cells with inferCNV mean absolute log2 burden at or above the global 75th percentile",
  cnv_burden_threshold_q75 = cnv_q75,
  candidate_T_cells_with_finite_CNV = nrow(candidate_cnv),
  selected_potential_malignant_T_cells = nrow(potential_cells),
  selected_percent = 100 * nrow(potential_cells) / nrow(candidate_cnv),
  stringsAsFactors = FALSE
)
utils::write.csv(
  selection_summary,
  file.path(final_csv, "potential_malignant_T_cell_definition.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

counts <- utils::read.csv(count_file, check.names = FALSE)
broad <- counts[counts$compartment == "broad_T_lineage_candidate_malignant", c("patient_pair", "sample_id", "condition_simple", "n_cells")]
names(broad)[4] <- "candidate_T_cells"
top25 <- counts[counts$compartment == "CNV_high_top25_candidate", c("patient_pair", "sample_id", "condition_simple", "n_cells")]
names(top25)[4] <- "potential_malignant_T_cells"
count_summary <- merge(broad, top25, by = c("patient_pair", "sample_id", "condition_simple"), all = FALSE)
count_summary$potential_malignant_T_cell_percent <- 100 * count_summary$potential_malignant_T_cells / count_summary$candidate_T_cells
count_summary <- count_summary[order(count_summary$patient_pair, count_summary$condition_simple), ]
utils::write.csv(
  count_summary,
  file.path(final_csv, "potential_malignant_T_cell_counts_and_proportions_by_sample.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

theme_clean <- function(base_size = 11) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(color = "black"),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    )
}

gsea_plot <- gsea
gsea_plot$pathway_label <- gsub("STAT5_", "", gsea_plot$pathway)
gsea_plot$pathway_label <- gsub("_", " ", gsea_plot$pathway_label)
gsea_plot$pathway_label <- factor(gsea_plot$pathway_label, levels = gsea_plot$pathway_label[order(gsea_plot$NES)])
gsea_plot$fdr_label <- sprintf("NES %.2f | FDR %.3f", gsea_plot$NES, gsea_plot$padj)
p_gsea <- ggplot2::ggplot(gsea_plot, ggplot2::aes(x = pathway_label, y = NES, fill = direction)) +
  ggplot2::geom_col(width = 0.7) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  ggplot2::geom_text(ggplot2::aes(label = fdr_label), hjust = -0.04, size = 3.5) +
  ggplot2::coord_flip(clip = "off") +
  ggplot2::scale_fill_manual(values = c(relapse_up = "#B2433F", relapse_down = "#3C6E71", flat = "grey70")) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.04, 0.42))) +
  ggplot2::labs(
    title = "STAT5 downstream gene sets in potential malignant T cells",
    subtitle = "CNV-high top 25%; relapse vs diagnosis; paired pseudobulk (N = 3)",
    x = NULL,
    y = "Normalized enrichment score (NES)"
  ) +
  theme_clean(11)
ggplot2::ggsave(
  file.path(final_png, "STAT5_downstream_gene_set_GSEA.png"),
  p_gsea,
  width = 8.2,
  height = 4.2,
  dpi = 300,
  bg = "white"
)

focus_names <- c(
  "Oxidative phosphorylation", "Mitophagy - animal", "Proteasome",
  "Chemical carcinogenesis - reactive oxygen species", "JAK-STAT signaling pathway",
  "Autophagy - animal"
)
kplot <- kegg_focus_top25[kegg_focus_top25$pathway_name %in% focus_names & kegg_focus_top25$tested, , drop = FALSE]
kplot$pathway_label <- factor(kplot$pathway_name, levels = kplot$pathway_name[order(kplot$NES)])
kplot$label <- sprintf("NES %.2f | FDR %.3f", kplot$NES, kplot$padj)
p_kegg <- ggplot2::ggplot(kplot, ggplot2::aes(x = pathway_label, y = NES, fill = direction)) +
  ggplot2::geom_col(width = 0.7) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  ggplot2::geom_text(ggplot2::aes(label = label), hjust = -0.04, size = 3.2) +
  ggplot2::coord_flip(clip = "off") +
  ggplot2::scale_fill_manual(values = c(relapse_up = "#B2433F", relapse_down = "#3C6E71", flat = "grey70")) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.35, 0.40))) +
  ggplot2::labs(
    title = "Key pathways in potential malignant T cells",
    subtitle = "CNV-high top 25%; relapse vs diagnosis; paired pseudobulk (N = 3)",
    x = NULL,
    y = "Normalized enrichment score (NES)"
  ) +
  theme_clean(10)
ggplot2::ggsave(
  file.path(final_png, "KEGG_GSEA_key_pathways.png"),
  p_kegg,
  width = 8.6,
  height = 5.2,
  dpi = 300,
  bg = "white"
)

custom_plot <- custom_mito_top25[custom_mito_top25$tested, , drop = FALSE]
custom_plot$pathway_label <- gsub("_", " ", custom_plot$pathway_name)
custom_plot$pathway_label <- factor(custom_plot$pathway_label, levels = custom_plot$pathway_label[order(custom_plot$NES)])
custom_plot$label <- sprintf("NES %.2f | FDR %.3f", custom_plot$NES, custom_plot$padj)
p_custom <- ggplot2::ggplot(custom_plot, ggplot2::aes(x = pathway_label, y = NES, fill = direction)) +
  ggplot2::geom_col(width = 0.7) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  ggplot2::geom_text(ggplot2::aes(label = label), hjust = -0.04, size = 3.2) +
  ggplot2::coord_flip(clip = "off") +
  ggplot2::scale_fill_manual(values = c(relapse_up = "#B2433F", relapse_down = "#3C6E71", flat = "grey70")) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.04, 0.43))) +
  ggplot2::labs(
    title = "Mitochondrial programs in potential malignant T cells",
    subtitle = "CNV-high top 25%; relapse vs diagnosis; paired pseudobulk (N = 3)",
    x = NULL,
    y = "Normalized enrichment score (NES)"
  ) +
  theme_clean(10)
ggplot2::ggsave(
  file.path(final_png, "custom_mitochondrial_gene_set_GSEA.png"),
  p_custom,
  width = 8.2,
  height = 5.0,
  dpi = 300,
  bg = "white"
)

stat5_ab$gene <- factor(stat5_ab$gene, levels = c("STAT5B", "STAT5A"))
stat5_ab$label <- sprintf("logFC %.2f\nnominal P %.3g", stat5_ab$logFC, stat5_ab$PValue)
p_ab <- ggplot2::ggplot(stat5_ab, ggplot2::aes(x = gene, y = logFC, fill = logFC > 0)) +
  ggplot2::geom_col(width = 0.62) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  ggplot2::geom_text(ggplot2::aes(label = label), vjust = ifelse(stat5_ab$logFC > 0, -0.25, 1.15), size = 3.5) +
  ggplot2::scale_fill_manual(values = c(`TRUE` = "#B2433F", `FALSE` = "#3C6E71")) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.18, 0.25))) +
  ggplot2::labs(
    title = "STAT5A and STAT5B expression in potential malignant T cells",
    subtitle = "CNV-high top 25%; relapse vs diagnosis; paired pseudobulk (N = 3)",
    x = NULL,
    y = "edgeR log2 fold change"
  ) +
  theme_clean(11)
ggplot2::ggsave(
  file.path(final_png, "STAT5A_STAT5B_pseudobulk_logFC.png"),
  p_ab,
  width = 6.2,
  height = 4.6,
  dpi = 300,
  bg = "white"
)

count_plot <- count_summary
count_plot$sample_label <- paste0(count_plot$patient_pair, " ", count_plot$sample_id)
count_plot$condition_simple <- factor(count_plot$condition_simple, levels = c("diagnosis", "relapse"))
count_plot$label_y <- count_plot$potential_malignant_T_cell_percent
count_plot$label_y[count_plot$condition_simple == "relapse" & count_plot$patient_pair == "P1"] <- 30.6
count_plot$label_y[count_plot$condition_simple == "relapse" & count_plot$patient_pair == "P3"] <- 35.0
p_prop <- ggplot2::ggplot(
  count_plot,
  ggplot2::aes(x = condition_simple, y = potential_malignant_T_cell_percent, group = patient_pair)
) +
  ggplot2::geom_line(color = "grey55", linewidth = 0.7) +
  ggplot2::geom_point(ggplot2::aes(shape = condition_simple), size = 3.2, color = "#B2433F", fill = "white", stroke = 1.1) +
  ggplot2::geom_text(
    ggplot2::aes(y = label_y, label = sprintf("%s: %.1f%%", patient_pair, potential_malignant_T_cell_percent)),
    hjust = ifelse(count_plot$condition_simple == "diagnosis", 1.12, -0.12),
    size = 3.4
  ) +
  ggplot2::scale_shape_manual(values = c(diagnosis = 21, relapse = 16)) +
  ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0.48, 0.48))) +
  ggplot2::labs(
    title = "Potential malignant T-cell proportion by patient",
    subtitle = "CNV-high top 25% among candidate malignant-like T-lineage cells",
    x = NULL,
    y = "Potential malignant T cells (%)"
  ) +
  theme_clean(11)
ggplot2::ggsave(
  file.path(final_png, "potential_malignant_T_cell_proportion_paired.png"),
  p_prop,
  width = 7.0,
  height = 5.2,
  dpi = 300,
  bg = "white"
)

cat("STAT5 GSEA complete.\n")
print(gsea[, c("pathway", "NES", "pval", "padj", "direction", "significant_FDR_0.05")])
