suppressPackageStartupMessages({
  library(DESeq2)
  library(fgsea)
  library(data.table)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

args <- commandArgs(trailingOnly = TRUE)
analysis_dir <- if (length(args) >= 1) args[[1]] else stop("Usage: Rscript analysis_gse218858_R.R <analysis_dir>")
analysis_dir <- normalizePath(analysis_dir, winslash = "/", mustWork = FALSE)

set.seed(218858)

download_dir <- file.path(analysis_dir, "downloads")
reference_dir <- file.path(analysis_dir, "reference")
processed_dir <- file.path(analysis_dir, "processed")
results_dir <- file.path(analysis_dir, "results")
de_dir <- file.path(results_dir, "de_results")
gsea_dir <- file.path(results_dir, "gsea")
plots_dir <- file.path(analysis_dir, "plots")
code_dir <- file.path(analysis_dir, "code_availability")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(code_dir, recursive = TRUE, showWarnings = FALSE)

main_contrast <- list(
  name = "R_DESeq2_primary_R2S5b_DN_vs_R2b_DN",
  case = "R2S5b_DN",
  control = "R2b_DN",
  label = "STAT5B N642H; Rag2-/- DN vs Rag2-/- DN"
)

key_genes <- c(
  "Stat5b", "Pim1", "Bcl2l1", "Socs2", "Cish", "Il2ra", "Il7r",
  "Lck", "Zap70", "Fyn", "Cd3d", "Cd3e", "Cd247", "Grap2",
  "Notch1", "Myc", "Mki67", "Dntt"
)

priority_terms <- c(
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_APOPTOSIS",
  "JAK-STAT signaling pathway",
  "T cell receptor signaling pathway",
  "Cytokine-cytokine receptor interaction",
  "Apoptosis",
  "Cell cycle",
  "Acute myeloid leukemia",
  "Human T-cell leukemia virus 1 infection"
)

strip_quotes <- function(x) {
  x <- trimws(x)
  sub('^"(.*)"$', "\\1", x)
}

parse_series_matrix <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE, encoding = "UTF-8")
  sample_lines <- lines[grepl("^!Sample_", lines)]
  rows <- lapply(sample_lines, function(line) strsplit(line, "\t", fixed = TRUE)[[1]])
  keys <- vapply(rows, `[[`, character(1), 1)
  vals <- lapply(rows, function(x) strip_quotes(x[-1]))
  sample_ids <- vals[[which(keys == "!Sample_geo_accession")[1]]]
  meta <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)

  for (i in seq_along(rows)) {
    key <- keys[[i]]
    value <- vals[[i]]
    if (length(value) != length(sample_ids)) next
    if (key == "!Sample_characteristics_ch1") {
      for (j in seq_along(value)) {
        item <- value[[j]]
        if (grepl(":", item, fixed = TRUE)) {
          parts <- strsplit(item, ":", fixed = TRUE)[[1]]
          k <- trimws(parts[[1]])
          v <- trimws(paste(parts[-1], collapse = ":"))
          k <- tolower(k)
          if (startsWith(k, "stage ")) k <- "stage"
          col <- make.names(k)
          if (!col %in% names(meta)) meta[[col]] <- NA_character_
          meta[[col]][j] <- v
        }
      }
    } else {
      clean_key <- tolower(sub("^!Sample_", "", key))
      meta[[make.names(clean_key)]] <- value
    }
  }

  meta$genotype_raw <- if ("genotype" %in% names(meta)) meta$genotype else ""
  meta$condition <- if ("condition" %in% names(meta)) meta$condition else ""
  meta$stage <- if ("stage" %in% names(meta)) meta$stage else sub("^.*_", "", meta$condition)
  meta$stage <- ifelse(meta$stage == "CD8", "SP8", meta$stage)
  meta$replicate_id <- if ("replicate.id" %in% names(meta)) meta$replicate.id else ""
  meta$name_stage_replicate <- if ("name.stage_replicate" %in% names(meta)) meta$name.stage_replicate else ""
  meta$rag2_status <- ifelse(grepl("Rag2", meta$genotype_raw, ignore.case = TRUE), "Rag2-/-", "Rag2+/+")
  meta$stat5b_n642h <- grepl("STAT5BN642H", meta$genotype_raw, ignore.case = TRUE)
  meta$stat5a_s710f <- grepl("Stat5aS710F", meta$genotype_raw, ignore.case = TRUE)
  meta$stat5_variant <- ifelse(meta$stat5b_n642h, "STAT5B_N642H", ifelse(meta$stat5a_s710f, "STAT5A_S710F", "none"))
  meta$genotype_clean <- dplyr::case_when(
    grepl("^R2S5b", meta$condition) ~ "STAT5B_N642H;Rag2-/-",
    grepl("^S5b", meta$condition) ~ "STAT5B_N642H",
    grepl("^R2S5a", meta$condition) ~ "STAT5A_S710F;Rag2-/-",
    grepl("^S5a", meta$condition) ~ "STAT5A_S710F",
    grepl("^R2", meta$condition) ~ "Rag2-/-",
    grepl("^WT", meta$condition) ~ "WT",
    TRUE ~ meta$genotype_raw
  )
  meta$is_main_focus <- meta$condition %in% c("R2b_DN", "R2S5b_DN")
  meta$is_context_focus <- meta$condition %in% c("R2b_DN", "R2S5b_DN", "R2S5b_DP", "R2S5b_CD8")
  meta$is_secondary_context <- meta$condition %in% c("WT_DN", "WT_DP", "WT_CD8", "S5b_DN", "S5b_DP", "S5b_CD8")
  meta$is_focus <- meta$is_context_focus | meta$is_secondary_context
  meta$sample_label <- paste0(meta$condition, "_", meta$replicate_id)
  meta <- meta[order(meta$sample_id), , drop = FALSE]
  meta
}

load_counts <- function(raw_tar, meta) {
  out_file <- file.path(processed_dir, "R_raw_count_matrix.csv.gz")
  if (file.exists(out_file)) {
    mat <- data.table::fread(out_file)
    ids <- mat[[1]]
    mat <- as.data.frame(mat[, -1])
    rownames(mat) <- ids
    return(as.matrix(mat))
  }

  extract_dir <- file.path(processed_dir, "raw_count_files")
  dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
  if (length(list.files(extract_dir, pattern = "\\.tsv\\.gz$")) == 0) {
    utils::untar(raw_tar, exdir = extract_dir)
  }

  files <- list.files(extract_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)
  names(files) <- sub("_.*$", "", basename(files))
  missing <- setdiff(meta$sample_id, names(files))
  if (length(missing) > 0) stop("Missing count files for: ", paste(missing, collapse = ", "))

  count_list <- lapply(meta$sample_id, function(sid) {
    dt <- data.table::fread(files[[sid]], header = FALSE, col.names = c("ensembl_id", sid))
    dt
  })
  merged <- Reduce(function(x, y) merge(x, y, by = "ensembl_id", all = TRUE), count_list)
  merged[is.na(merged)] <- 0
  data.table::fwrite(merged, out_file)
  mat <- as.data.frame(merged[, -1])
  rownames(mat) <- merged$ensembl_id
  as.matrix(mat)
}

annotate_genes <- function(ensembl_ids) {
  ids <- sub("\\..*$", "", ensembl_ids)
  suppressMessages({
    ann <- AnnotationDbi::select(
      org.Mm.eg.db,
      keys = unique(ids),
      keytype = "ENSEMBL",
      columns = c("SYMBOL", "GENENAME", "ENTREZID")
    )
  })
  ann <- ann[!duplicated(ann$ENSEMBL), , drop = FALSE]
  ann <- data.frame(ensembl_id = ids, stringsAsFactors = FALSE) |>
    left_join(ann, by = c("ensembl_id" = "ENSEMBL")) |>
    rename(gene_symbol = SYMBOL, gene_name = GENENAME, entrez_id = ENTREZID)
  ann
}

log2_cpm <- function(counts) {
  lib <- colSums(counts)
  cpm <- sweep(counts, 2, lib, "/") * 1e6
  log2(cpm + 1)
}

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  split_lines <- strsplit(lines, "\t", fixed = TRUE)
  pathways <- lapply(split_lines, function(x) unique(x[-c(1, 2)]))
  names(pathways) <- vapply(split_lines, `[[`, character(1), 1)
  pathways
}

map_pathways_to_rank_symbols <- function(pathways, rank_symbols) {
  all_genes <- unique(unlist(pathways, use.names = FALSE))
  upper_map <- rank_symbols[!duplicated(toupper(rank_symbols))]
  names(upper_map) <- toupper(upper_map)

  direct_map <- rank_symbols
  names(direct_map) <- rank_symbols
  case_map <- upper_map[toupper(all_genes)]
  names(case_map) <- all_genes

  mapped <- direct_map[all_genes]
  names(mapped) <- all_genes
  missing <- is.na(mapped)
  mapped[missing] <- case_map[names(mapped)[missing]]

  still_missing <- names(mapped)[is.na(mapped)]
  if (length(still_missing) > 0) {
    alias <- tryCatch(
      suppressMessages(
        AnnotationDbi::select(
          org.Mm.eg.db,
          keys = unique(still_missing),
          keytype = "ALIAS",
          columns = c("SYMBOL")
        )
      ),
      error = function(e) data.frame(ALIAS = character(), SYMBOL = character())
    )
    alias <- alias[!is.na(alias$SYMBOL) & alias$SYMBOL %in% rank_symbols, , drop = FALSE]
    alias <- alias[!duplicated(alias$ALIAS), , drop = FALSE]
    alias_map <- alias$SYMBOL
    names(alias_map) <- alias$ALIAS
    idx <- names(mapped) %in% names(alias_map) & is.na(mapped)
    mapped[idx] <- alias_map[names(mapped)[idx]]
  }

  lapply(pathways, function(genes) {
    out <- unname(mapped[genes])
    unique(out[!is.na(out) & out != ""])
  })
}

collapse_leading_edge <- function(x) {
  vapply(x, function(y) paste(y, collapse = ";"), character(1))
}

run_fgsea_collection <- function(rank_stats, pathways, collection, permutations = 10000) {
  set.seed(218858)
  res <- suppressWarnings(
    fgsea::fgsea(
      pathways = pathways,
      stats = rank_stats,
      minSize = 5,
      maxSize = 2000,
      nperm = permutations
    )
  )
  res <- as.data.frame(res)
  if (!"log2err" %in% names(res)) res$log2err <- NA_real_
  if ("leadingEdge" %in% names(res)) res$leadingEdge <- collapse_leading_edge(res$leadingEdge)
  res$collection <- collection
  res$permutation_num <- permutations
  res <- res |>
    dplyr::select(collection, pathway, pval, padj, log2err, ES, NES, size, leadingEdge, permutation_num) |>
    arrange(padj, pval)
  res
}

save_plot <- function(filename, plot, width = 8, height = 6) {
  ggsave(file.path(plots_dir, filename), plot, width = width, height = height, dpi = 240, bg = "white")
}

message("Parsing metadata...")
series_matrix <- file.path(download_dir, "GSE218858_series_matrix.txt.gz")
raw_tar <- file.path(download_dir, "GSE218858_RAW.tar")
meta <- parse_series_matrix(series_matrix)
data.table::fwrite(meta, file.path(processed_dir, "R_sample_annotation.csv"))

message("Loading raw counts...")
counts <- load_counts(raw_tar, meta)
storage.mode(counts) <- "integer"
counts <- counts[, meta$sample_id, drop = FALSE]

message("Annotating genes...")
annotation <- annotate_genes(rownames(counts))
data.table::fwrite(annotation, file.path(processed_dir, "R_gene_annotation_orgMm.csv"))

message("QC and normalized matrices...")
lcpm <- log2_cpm(counts)
data.table::fwrite(data.frame(ensembl_id = rownames(lcpm), lcpm, check.names = FALSE), file.path(processed_dir, "R_log2_cpm_matrix.csv.gz"))
qc <- data.frame(
  sample_id = colnames(counts),
  library_size = colSums(counts),
  detected_genes_count_gt0 = colSums(counts > 0),
  genes_count_ge10 = colSums(counts >= 10),
  stringsAsFactors = FALSE
) |>
  left_join(meta, by = "sample_id")
data.table::fwrite(qc, file.path(processed_dir, "R_sample_qc_metrics.csv"))

qc_plot <- ggplot(qc |> filter(is_focus), aes(x = sample_label, y = library_size, fill = genotype_clean)) +
  geom_col() +
  facet_grid(. ~ stage, scales = "free_x", space = "free_x") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "bottom") +
  labs(title = "Library size - focus/context samples", x = NULL, y = "Total assigned counts", fill = "Genotype")
save_plot("R_sample_qc_library_size.png", qc_plot, width = 10, height = 5)

message("Running DESeq2 main contrast...")
main_meta <- meta |> filter(condition %in% c(main_contrast$control, main_contrast$case))
main_counts <- counts[, main_meta$sample_id, drop = FALSE]
keep <- rowSums(main_counts) >= 10 & rowSums(main_counts >= 5) >= 2
main_counts <- main_counts[keep, , drop = FALSE]
rownames(main_meta) <- main_meta$sample_id
main_meta$condition <- relevel(factor(main_meta$condition), ref = main_contrast$control)
dds <- DESeqDataSetFromMatrix(
  countData = round(main_counts),
  colData = main_meta,
  design = ~ condition
)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("condition", main_contrast$case, main_contrast$control), alpha = 0.05)
res_df <- as.data.frame(res) |>
  rownames_to_column("ensembl_id") |>
  left_join(annotation, by = "ensembl_id") |>
  mutate(
    contrast = main_contrast$name,
    contrast_label = main_contrast$label,
    case_condition = main_contrast$case,
    control_condition = main_contrast$control,
    method = "R DESeq2 negative-binomial Wald test",
    significant_fdr_0_05 = !is.na(padj) & padj < 0.05,
    significant_fdr_0_05_abs_lfc_1 = significant_fdr_0_05 & abs(log2FoldChange) >= 1,
    rank_metric = ifelse(is.na(stat), sign(log2FoldChange) * -log10(pmax(pvalue, 1e-300)), stat)
  ) |>
  arrange(padj, pvalue)
data.table::fwrite(res_df, file.path(de_dir, paste0(main_contrast$name, ".csv")))

ranked <- res_df |>
  filter(!is.na(gene_symbol), gene_symbol != "", !is.na(rank_metric)) |>
  group_by(gene_symbol) |>
  slice_max(order_by = abs(rank_metric), n = 1, with_ties = FALSE) |>
  ungroup() |>
  arrange(desc(rank_metric))
data.table::fwrite(ranked, file.path(results_dir, "R_ranked_genes_main_contrast.csv"))
data.table::fwrite(ranked |> dplyr::select(gene_symbol, rank_metric), file.path(results_dir, "R_ranked_genes_main_contrast.rnk"), sep = "\t", col.names = FALSE)

de_summary <- data.frame(
  contrast = main_contrast$name,
  contrast_label = main_contrast$label,
  method = "R DESeq2 negative-binomial Wald test",
  genes_tested = nrow(res_df),
  fdr_0_05 = sum(res_df$significant_fdr_0_05, na.rm = TRUE),
  fdr_0_05_up = sum(res_df$significant_fdr_0_05 & res_df$log2FoldChange > 0, na.rm = TRUE),
  fdr_0_05_down = sum(res_df$significant_fdr_0_05 & res_df$log2FoldChange < 0, na.rm = TRUE),
  fdr_0_05_abs_lfc_1 = sum(res_df$significant_fdr_0_05_abs_lfc_1, na.rm = TRUE),
  stringsAsFactors = FALSE
)
data.table::fwrite(de_summary, file.path(results_dir, "R_DESeq2_main_contrast_summary.csv"))

message("PCA and heatmaps...")
focus_meta <- meta |> filter(is_focus)
focus_counts <- counts[, focus_meta$sample_id, drop = FALSE]
focus_keep <- rowSums(focus_counts) >= 10
focus_dds <- DESeqDataSetFromMatrix(
  countData = round(focus_counts[focus_keep, , drop = FALSE]),
  colData = as.data.frame(focus_meta |> column_to_rownames("sample_id")),
  design = ~ condition
)
focus_vst <- vst(focus_dds, blind = TRUE)
vst_mat <- assay(focus_vst)
top_var <- head(order(matrixStats::rowVars(vst_mat), decreasing = TRUE), min(5000, nrow(vst_mat)))
pca <- prcomp(t(vst_mat[top_var, , drop = FALSE]), scale. = FALSE)
pca_df <- as.data.frame(pca$x[, 1:3, drop = FALSE]) |>
  rownames_to_column("sample_id") |>
  left_join(meta, by = "sample_id")
data.table::fwrite(pca_df, file.path(processed_dir, "R_PCA_coordinates_focus_samples.csv"))
pct <- round((pca$sdev^2 / sum(pca$sdev^2))[1:2] * 100, 1)
pca_plot <- ggplot(pca_df, aes(PC1, PC2, color = genotype_clean, shape = stage, label = sample_label)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 2.4, max.overlaps = 50, show.legend = FALSE) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom") +
  labs(title = "R/DESeq2 VST PCA - focus/context samples", x = paste0("PC1 (", pct[[1]], "%)"), y = paste0("PC2 (", pct[[2]], "%)"), color = "Genotype")
save_plot("R_PCA_focus_samples.png", pca_plot, width = 9, height = 6.5)

volcano <- res_df |>
  mutate(
    neg_log10_padj = -log10(pmax(padj, 1e-300)),
    category = case_when(
      !is.na(padj) & padj < 0.05 & log2FoldChange >= 1 ~ "up",
      !is.na(padj) & padj < 0.05 & log2FoldChange <= -1 ~ "down",
      TRUE ~ "not significant"
    ),
    label_gene = ifelse(gene_symbol %in% key_genes & !is.na(padj) & padj < 0.1, gene_symbol, NA)
  )
volcano_plot <- ggplot(volcano, aes(log2FoldChange, neg_log10_padj, color = category)) +
  geom_point(size = 0.7, alpha = 0.7) +
  ggrepel::geom_text_repel(aes(label = label_gene), size = 2.5, max.overlaps = 50, na.rm = TRUE, show.legend = FALSE) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("up" = "#B23A48", "down" = "#2F6F9F", "not significant" = "#B8B8B8")) +
  theme_bw(base_size = 10) +
  labs(title = main_contrast$label, x = "log2 fold change", y = "-log10(FDR)", color = NULL)
save_plot("R_volcano_R2S5b_DN_vs_R2b_DN.png", volcano_plot, width = 7.5, height = 5.8)

top_de_ids <- res_df |>
  filter(!is.na(padj)) |>
  arrange(padj) |>
  slice_head(n = 80) |>
  pull(ensembl_id)
heat_mat <- vst_mat[intersect(top_de_ids, rownames(vst_mat)), focus_meta$sample_id, drop = FALSE]
heat_mat <- t(scale(t(heat_mat)))
row_labels <- annotation$gene_symbol[match(rownames(heat_mat), annotation$ensembl_id)]
row_labels[is.na(row_labels) | row_labels == ""] <- rownames(heat_mat)[is.na(row_labels) | row_labels == ""]
rownames(heat_mat) <- make.unique(row_labels)
ann_col <- focus_meta |>
  dplyr::select(sample_id, condition, genotype_clean, stage) |>
  as.data.frame()
rownames(ann_col) <- ann_col$sample_id
ann_col$sample_id <- NULL
png(file.path(plots_dir, "R_top80_DE_genes_heatmap_focus_samples.png"), width = 2100, height = 2400, res = 240)
pheatmap(
  heat_mat,
  annotation_col = ann_col,
  show_colnames = TRUE,
  fontsize_row = 5.5,
  fontsize_col = 7,
  main = "Top DE genes, R2S5b_DN vs R2b_DN"
)
dev.off()

best_key <- annotation |>
  filter(gene_symbol %in% key_genes) |>
  left_join(res_df |> dplyr::select(ensembl_id, baseMean), by = "ensembl_id") |>
  arrange(gene_symbol, desc(baseMean)) |>
  group_by(gene_symbol) |>
  slice_head(n = 1) |>
  ungroup()
key_expr <- vst_mat[intersect(best_key$ensembl_id, rownames(vst_mat)), focus_meta$sample_id, drop = FALSE]
key_symbols <- best_key$gene_symbol[match(rownames(key_expr), best_key$ensembl_id)]
rownames(key_expr) <- key_symbols
key_z <- t(scale(t(key_expr)))
png(file.path(plots_dir, "R_selected_STAT5_TCR_genes_heatmap.png"), width = 1900, height = 1200, res = 240)
pheatmap(
  key_z,
  annotation_col = ann_col,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 7,
  main = "Selected STAT5/TCR/proliferation/survival genes"
)
dev.off()

key_de <- res_df |> filter(gene_symbol %in% key_genes) |> arrange(padj, gene_symbol)
data.table::fwrite(key_de, file.path(results_dir, "R_selected_STAT5_TCR_genes_DESeq2.csv"))

message("Running fgsea Hallmark/KEGG/GO BP...")
rank_stats <- ranked$rank_metric
names(rank_stats) <- ranked$gene_symbol
rank_stats <- sort(rank_stats, decreasing = TRUE)

hallmark <- read_gmt(file.path(reference_dir, "mh.all.v2024.1.Mm.symbols.gmt"))
kegg <- read_gmt(file.path(reference_dir, "KEGG_2019_Mouse.enrichr.gmt"))
kegg <- map_pathways_to_rank_symbols(kegg, names(rank_stats))
data.table::fwrite(
  data.frame(
    pathway = names(kegg),
    mapped_gene_count = vapply(kegg, length, integer(1))
  ),
  file.path(gsea_dir, "R_KEGG_2019_Mouse_symbol_mapping_counts.csv")
)
gobp <- read_gmt(file.path(reference_dir, "m5.go.bp.v2024.1.Mm.symbols.gmt"))

gsea_h <- run_fgsea_collection(rank_stats, hallmark, "Hallmark_MSigDB_2024_1_Mm", permutations = 10000)
gsea_k <- run_fgsea_collection(rank_stats, kegg, "KEGG_2019_Mouse_Enrichr_mapped", permutations = 10000)
gsea_go <- run_fgsea_collection(rank_stats, gobp, "GO_BP_MSigDB_2024_1_Mm", permutations = 10000)
gsea_all <- bind_rows(gsea_h, gsea_k, gsea_go) |>
  mutate(
    contrast = main_contrast$name,
    contrast_label = main_contrast$label,
    significant_fdr_0_05 = !is.na(padj) & padj < 0.05,
    nominal_fdr_0_25 = !is.na(padj) & padj < 0.25
  ) |>
  dplyr::select(contrast, contrast_label, collection, pathway, pval, padj, log2err, ES, NES, size, leadingEdge, permutation_num, significant_fdr_0_05, nominal_fdr_0_25) |>
  arrange(collection, padj, pval)

data.table::fwrite(gsea_all, file.path(gsea_dir, "R_fgsea_all_pathways_Hallmark_KEGG_GO_BP.csv"))
data.table::fwrite(gsea_all |> filter(significant_fdr_0_05), file.path(gsea_dir, "R_fgsea_significant_pathways_fdr_0_05.csv"))
data.table::fwrite(gsea_all |> filter(nominal_fdr_0_25), file.path(gsea_dir, "R_fgsea_nominal_pathways_fdr_0_25.csv"))
for (coll in unique(gsea_all$collection)) {
  short <- gsub("[^A-Za-z0-9_]+", "_", coll)
  sub <- gsea_all |> filter(.data$collection == coll)
  data.table::fwrite(sub, file.path(gsea_dir, paste0("R_fgsea_", short, "_all_pathways.csv")))
  data.table::fwrite(sub |> filter(significant_fdr_0_05), file.path(gsea_dir, paste0("R_fgsea_", short, "_significant_pathways_fdr_0_05.csv")))
}

priority <- gsea_all |> filter(pathway %in% priority_terms) |> arrange(padj, desc(abs(NES)))
data.table::fwrite(priority, file.path(gsea_dir, "R_fgsea_priority_terms_main_contrast.csv"))

priority_plot_df <- priority |> mutate(neg_log10_fdr = -log10(pmax(padj, 1e-300)))
priority_plot <- ggplot(priority_plot_df, aes(NES, reorder(pathway, NES), size = neg_log10_fdr, color = collection)) +
  geom_point(alpha = 0.85) +
  geom_vline(xintercept = 0, color = "grey50") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom") +
  labs(title = "Priority pathway fgsea - R main contrast", x = "NES", y = NULL, size = "-log10(FDR)")
save_plot("R_fgsea_priority_pathways_dotplot.png", priority_plot, width = 10, height = 6.5)

top_gsea <- gsea_all |>
  group_by(collection) |>
  arrange(padj) |>
  slice_head(n = 12) |>
  ungroup() |>
  mutate(pathway_short = substr(gsub("^HALLMARK_", "", pathway), 1, 60))
top_gsea_plot <- ggplot(top_gsea, aes(NES, reorder(pathway_short, NES), fill = collection)) +
  geom_col() +
  facet_wrap(~ collection, scales = "free_y") +
  geom_vline(xintercept = 0, color = "grey50") +
  theme_bw(base_size = 9) +
  theme(legend.position = "none") +
  labs(title = "Top fgsea pathways - R main contrast", x = "NES", y = NULL)
save_plot("R_fgsea_top_pathways_barplot.png", top_gsea_plot, width = 12, height = 8)

message("Pathway scores...")
pathway_score_terms <- priority_terms[priority_terms %in% names(c(hallmark, kegg))]
symbol_vst <- data.frame(ensembl_id = rownames(vst_mat), vst_mat, check.names = FALSE) |>
  left_join(annotation, by = "ensembl_id") |>
  filter(!is.na(gene_symbol), gene_symbol != "") |>
  mutate(mean_expr = rowMeans(across(all_of(focus_meta$sample_id)))) |>
  arrange(gene_symbol, desc(mean_expr)) |>
  group_by(gene_symbol) |>
  slice_head(n = 1) |>
  ungroup()
symbol_mat <- as.matrix(symbol_vst[, focus_meta$sample_id, drop = FALSE])
rownames(symbol_mat) <- symbol_vst$gene_symbol
symbol_z <- t(scale(t(symbol_mat)))
all_sets <- c(hallmark, kegg)
score_rows <- list()
for (term in pathway_score_terms) {
  genes <- intersect(all_sets[[term]], rownames(symbol_z))
  if (length(genes) < 3) next
  vals <- colMeans(symbol_z[genes, , drop = FALSE], na.rm = TRUE)
  score_rows[[term]] <- data.frame(sample_id = names(vals), pathway = term, score = as.numeric(vals), matched_genes = length(genes))
}
scores <- bind_rows(score_rows) |> left_join(meta, by = "sample_id")
data.table::fwrite(scores, file.path(results_dir, "R_selected_pathway_scores.csv"))
score_plot <- ggplot(scores |> filter(is_focus), aes(condition, score, fill = stage)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.65) +
  geom_jitter(width = 0.12, size = 1.5) +
  facet_wrap(~ pathway, scales = "free_y", ncol = 3) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "bottom") +
  labs(title = "Selected pathway scores - focus/context samples", x = NULL, y = "Mean z-scored VST expression")
save_plot("R_selected_pathway_scores_boxplots.png", score_plot, width = 11, height = 8)

message("Comparing R and prior Python result if available...")
python_dir <- file.path(dirname(analysis_dir), "GSE218858_STAT5B_N642H_analysis")
py_de_file <- file.path(python_dir, "results", "de_results", "primary_R2S5b_DN_vs_R2b_DN.csv")
comparison_summary_file <- file.path(results_dir, "R_vs_Python_main_contrast_summary.csv")
comparison_flag_file <- file.path(results_dir, "R_vs_Python_similarity_passed.txt")
comparison_note <- "Prior Python result was not found; no comparison performed."
similar_to_python <- if (file.exists(comparison_flag_file)) {
  isTRUE(toupper(trimws(readLines(comparison_flag_file, warn = FALSE)[1])) == "TRUE")
} else {
  FALSE
}
if (file.exists(py_de_file)) {
  py <- data.table::fread(py_de_file)
  cmp <- res_df |>
    dplyr::select(ensembl_id, R_log2FoldChange = log2FoldChange, R_stat = stat, R_padj = padj, gene_symbol) |>
    inner_join(
      py |> dplyr::select(ensembl_id, Python_log2FoldChange = log2FoldChange, Python_stat = stat, Python_padj = padj),
      by = "ensembl_id"
    )
  lfc_cor <- suppressWarnings(cor(cmp$R_log2FoldChange, cmp$Python_log2FoldChange, method = "spearman", use = "complete.obs"))
  stat_cor <- suppressWarnings(cor(cmp$R_stat, cmp$Python_stat, method = "spearman", use = "complete.obs"))
  sign_concordance <- mean(sign(cmp$R_log2FoldChange) == sign(cmp$Python_log2FoldChange), na.rm = TRUE)
  r_sig <- cmp$ensembl_id[!is.na(cmp$R_padj) & cmp$R_padj < 0.05]
  py_sig <- cmp$ensembl_id[!is.na(cmp$Python_padj) & cmp$Python_padj < 0.05]
  jaccard_sig <- length(intersect(r_sig, py_sig)) / length(union(r_sig, py_sig))
  cmp_summary <- data.frame(
    shared_genes = nrow(cmp),
    spearman_log2fc = lfc_cor,
    spearman_stat = stat_cor,
    log2fc_sign_concordance = sign_concordance,
    significant_gene_jaccard_fdr_0_05 = jaccard_sig
  )
  similar_to_python <- isTRUE(lfc_cor >= 0.97 && stat_cor >= 0.95 && sign_concordance >= 0.95)
  cmp_summary$similar_to_python_threshold <- similar_to_python
  data.table::fwrite(cmp, file.path(results_dir, "R_vs_Python_main_contrast_gene_level_comparison.csv"))
  data.table::fwrite(cmp_summary, file.path(results_dir, "R_vs_Python_main_contrast_summary.csv"))
  comparison_note <- paste0(
    "R-vs-Python comparison for the main contrast: shared genes = ", nrow(cmp),
    ", Spearman log2FC = ", round(lfc_cor, 4),
    ", Spearman statistic = ", round(stat_cor, 4),
    ", log2FC sign concordance = ", round(sign_concordance, 4),
    ", FDR<0.05 Jaccard = ", round(jaccard_sig, 4),
    ". Similarity threshold passed = ", similar_to_python, "."
  )
} else if (file.exists(comparison_summary_file)) {
  cmp_summary <- data.table::fread(comparison_summary_file)
  comparison_note <- paste0(
    "Prior Python result folder was already removed, so the existing comparison summary was retained: shared genes = ",
    cmp_summary$shared_genes[[1]],
    ", Spearman log2FC = ", round(cmp_summary$spearman_log2fc[[1]], 4),
    ", Spearman statistic = ", round(cmp_summary$spearman_stat[[1]], 4),
    ", log2FC sign concordance = ", round(cmp_summary$log2fc_sign_concordance[[1]], 4),
    ", FDR<0.05 Jaccard = ", round(cmp_summary$significant_gene_jaccard_fdr_0_05[[1]], 4),
    ". Similarity threshold passed = ", similar_to_python, "."
  )
}

message("Writing report...")
top_up <- res_df |> filter(!is.na(padj)) |> arrange(desc(log2FoldChange)) |> slice_head(n = 15) |> pull(gene_symbol)
top_down <- res_df |> filter(!is.na(padj)) |> arrange(log2FoldChange) |> slice_head(n = 15) |> pull(gene_symbol)
priority_text <- priority |>
  dplyr::select(collection, pathway, NES, pval, padj) |>
  arrange(padj) |>
  mutate(across(c(NES, pval, padj), ~ signif(.x, 4)))
key_text <- key_de |>
  dplyr::select(gene_symbol, log2FoldChange, stat, pvalue, padj) |>
  arrange(padj) |>
  mutate(across(c(log2FoldChange, stat, pvalue, padj), ~ signif(.x, 4)))

report <- c(
  "# GSE218858 R/DESeq2 main contrast analysis",
  "",
  "## Scope",
  "",
  "This R version focuses on the primary mechanistic comparison: STAT5B N642H; Rag2-/- DN thymocytes (`R2S5b_DN`) versus Rag2-/- DN thymocytes (`R2b_DN`). DP/SP8 and WT/S5b samples are retained only for QC/PCA/context plots.",
  "",
  "This dataset should be interpreted as a mouse mechanistic immature T-lineage / T-ALL-like model. It does not prove bortezomib resistance, clinical relapse, or patient response biology.",
  "",
  "## Differential Expression",
  "",
  paste0("- Method: R DESeq2 ", as.character(packageVersion("DESeq2")), " negative-binomial Wald test"),
  paste0("- Genes tested: ", de_summary$genes_tested),
  paste0("- FDR < 0.05 genes: ", de_summary$fdr_0_05, " (up = ", de_summary$fdr_0_05_up, ", down = ", de_summary$fdr_0_05_down, ")"),
  paste0("- FDR < 0.05 and |log2FC| >= 1 genes: ", de_summary$fdr_0_05_abs_lfc_1),
  paste0("- Top up genes by log2FC: ", paste(na.omit(top_up), collapse = "; ")),
  paste0("- Top down genes by log2FC: ", paste(na.omit(top_down), collapse = "; ")),
  "",
  "## Priority Pathways",
  "",
  paste(capture.output(print(priority_text, row.names = FALSE)), collapse = "\n"),
  "",
  "## Selected STAT5/TCR/Survival/Proliferation Genes",
  "",
  paste(capture.output(print(key_text, row.names = FALSE)), collapse = "\n"),
  "",
  "## R vs Python",
  "",
  comparison_note,
  "",
  "## Main Output Files",
  "",
  "- `processed/R_sample_annotation.csv`",
  "- `processed/R_sample_qc_metrics.csv`",
  "- `results/de_results/R_DESeq2_primary_R2S5b_DN_vs_R2b_DN.csv`",
  "- `results/R_ranked_genes_main_contrast.csv`",
  "- `results/gsea/R_fgsea_all_pathways_Hallmark_KEGG_GO_BP.csv`",
  "- `results/gsea/R_fgsea_significant_pathways_fdr_0_05.csv`",
  "- `results/gsea/R_fgsea_priority_terms_main_contrast.csv`",
  "- `results/R_selected_STAT5_TCR_genes_DESeq2.csv`",
  "- `plots/*.png`",
  "",
  "## Session Info",
  "",
  "```",
  paste(capture.output(sessionInfo()), collapse = "\n"),
  "```"
)
writeLines(report, file.path(analysis_dir, "GSE218858_R_DESeq2_main_contrast_report.md"), useBytes = TRUE)
writeLines(as.character(similar_to_python), file.path(results_dir, "R_vs_Python_similarity_passed.txt"), useBytes = TRUE)

message("Done.")
