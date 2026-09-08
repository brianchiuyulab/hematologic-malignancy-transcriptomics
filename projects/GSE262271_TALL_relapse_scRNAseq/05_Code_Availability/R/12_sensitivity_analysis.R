options(stringsAsFactors = FALSE)
options(timeout = 300)

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
min_cells_per_sample <- as.integer(get_arg("--min-cells-per-sample", "20"))
if (is.na(min_cells_per_sample) || min_cells_per_sample < 1) min_cells_per_sample <- 20

paths <- list(
  objects = file.path(project, "processed_objects"),
  results = file.path(project, "results"),
  figures = file.path(project, "figures"),
  logs = file.path(project, "logs"),
  source = file.path(project, "source_metadata")
)
paths$sensitivity_results <- file.path(paths$results, "sensitivity_mitophagy_gsea")
paths$sensitivity_figures <- file.path(paths$figures, "sensitivity_mitophagy_gsea")
invisible(lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))

user_lib <- Sys.getenv("R_LIBS_USER", unset = "C:/Users/User/Documents/R/win-library/4.4")
if (dir.exists(user_lib)) .libPaths(unique(c(normalizePath(user_lib, winslash = "/"), .libPaths())))

load_or_stop <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
  suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))
}
load_or_stop(c(
  "Seurat", "SeuratObject", "Matrix", "dplyr", "tidyr", "tibble", "ggplot2",
  "edgeR", "fgsea", "AnnotationDbi", "org.Hs.eg.db", "GO.db"
))

log_file <- file.path(paths$logs, "sensitivity_mitophagy_gsea.log")
log_message <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

write_csv_file <- function(x, path) {
  utils::write.csv(as.data.frame(x), path, row.names = FALSE, fileEncoding = "UTF-8")
}

theme_pub <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(color = "black"),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "grey70"),
      legend.key = ggplot2::element_blank()
    )
}

slugify <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  substr(x, 1, 140)
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

download_if_missing <- function(url, dest) {
  if (!file.exists(dest) || file.info(dest)$size == 0) {
    log_message("Downloading", url)
    utils::download.file(url, dest, mode = "wb", quiet = TRUE)
  }
  dest
}

read_kegg_two_col <- function(path, col_names = c("V1", "V2")) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  split_once <- regexpr("\t", lines, fixed = TRUE)
  out <- data.frame(
    V1 = ifelse(split_once > 0, substr(lines, 1, split_once - 1), lines),
    V2 = ifelse(split_once > 0, substr(lines, split_once + 1, nchar(lines)), ""),
    stringsAsFactors = FALSE
  )
  stats::setNames(out, col_names)
}

build_kegg_gene_sets <- function(paths) {
  kegg_files <- list(
    pathways = download_if_missing("https://rest.kegg.jp/list/pathway/hsa", file.path(paths$source, "kegg_hsa_pathway_list.tsv")),
    links = download_if_missing("https://rest.kegg.jp/link/pathway/hsa", file.path(paths$source, "kegg_hsa_pathway_gene_links.tsv")),
    genes = download_if_missing("https://rest.kegg.jp/list/hsa", file.path(paths$source, "kegg_hsa_gene_list.tsv"))
  )

  pathway_tbl <- read_kegg_two_col(kegg_files$pathways, c("pathway_id", "pathway_name"))
  pathway_tbl$pathway_id <- sub("^path:", "", pathway_tbl$pathway_id)
  pathway_tbl$pathway_name <- sub(" - Homo sapiens \\(human\\)$", "", pathway_tbl$pathway_name)

  gene_tbl <- read_kegg_two_col(kegg_files$genes, c("kegg_gene", "description"))
  gene_tbl$kegg_gene <- sub("^hsa:", "", gene_tbl$kegg_gene)
  gene_symbol_field <- ifelse(
    grepl("\t", gene_tbl$description, fixed = TRUE),
    sub("^[^\t]*\t[^\t]*\t", "", gene_tbl$description),
    gene_tbl$description
  )
  gene_tbl$primary_symbol <- trimws(sub(",.*$", "", sub(";.*$", "", gene_symbol_field)))
  gene_tbl <- gene_tbl[nzchar(gene_tbl$primary_symbol), c("kegg_gene", "primary_symbol")]
  gene_tbl <- gene_tbl[!duplicated(gene_tbl$kegg_gene), ]

  links_raw <- read_kegg_two_col(kegg_files$links, c("V1", "V2"))
  if (grepl("^path:", links_raw$V1[[1]])) {
    links <- data.frame(
      pathway_id = sub("^path:", "", links_raw$V1),
      kegg_gene = sub("^hsa:", "", links_raw$V2),
      stringsAsFactors = FALSE
    )
  } else {
    links <- data.frame(
      kegg_gene = sub("^hsa:", "", links_raw$V1),
      pathway_id = sub("^path:", "", links_raw$V2),
      stringsAsFactors = FALSE
    )
  }

  links_symbols <- dplyr::inner_join(links, gene_tbl, by = "kegg_gene", relationship = "many-to-one") |>
    dplyr::rename(symbol = primary_symbol) |>
    dplyr::left_join(pathway_tbl, by = "pathway_id") |>
    dplyr::mutate(pathway = paste0(pathway_id, " ", pathway_name))

  pathway_list <- links_symbols |>
    dplyr::filter(!is.na(pathway_name), nzchar(symbol)) |>
    dplyr::group_by(pathway) |>
    dplyr::summarise(genes = list(sort(unique(symbol))), .groups = "drop")

  audit <- pathway_tbl |>
    dplyr::mutate(
      gene_set_source = "KEGG",
      gene_set_id = pathway_id,
      pathway = paste0(pathway_id, " ", pathway_name)
    ) |>
    dplyr::left_join(
      links_symbols |>
        dplyr::group_by(pathway_id) |>
        dplyr::summarise(
          linked_genes_n = dplyr::n_distinct(symbol),
          linked_genes = paste(sort(unique(symbol)), collapse = ";"),
          .groups = "drop"
        ),
      by = "pathway_id"
    ) |>
    dplyr::mutate(
      linked_genes_n = ifelse(is.na(linked_genes_n), 0L, linked_genes_n),
      linked_genes = ifelse(is.na(linked_genes), "", linked_genes)
    )

  list(pathways = stats::setNames(pathway_list$genes, pathway_list$pathway), audit = audit)
}

build_go_bp_gene_sets <- function(symbol_universe) {
  log_message("Building GO BP gene sets from org.Hs.eg.db/GO.db")
  valid_symbols <- intersect(symbol_universe, AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "SYMBOL"))
  go_map <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = valid_symbols,
    columns = c("GOALL", "ONTOLOGYALL"),
    keytype = "SYMBOL"
  ))
  go_map <- go_map |>
    dplyr::filter(!is.na(GOALL), ONTOLOGYALL == "BP") |>
    dplyr::distinct(SYMBOL, GOALL)
  go_terms <- suppressMessages(AnnotationDbi::select(
    GO.db::GO.db,
    keys = unique(go_map$GOALL),
    columns = c("TERM", "ONTOLOGY"),
    keytype = "GOID"
  ))
  go_terms <- go_terms |>
    dplyr::filter(ONTOLOGY == "BP") |>
    dplyr::distinct(GOID, TERM)
  go_joined <- go_map |>
    dplyr::inner_join(go_terms, by = c("GOALL" = "GOID")) |>
    dplyr::mutate(pathway = paste0(GOALL, " ", TERM))

  go_list <- go_joined |>
    dplyr::group_by(pathway) |>
    dplyr::summarise(genes = list(sort(unique(SYMBOL))), .groups = "drop")

  audit <- go_joined |>
    dplyr::group_by(GOALL, TERM, pathway) |>
    dplyr::summarise(
      gene_set_source = "GO_BP",
      gene_set_id = dplyr::first(GOALL),
      pathway_name = dplyr::first(TERM),
      linked_genes_n = dplyr::n_distinct(SYMBOL),
      linked_genes = paste(sort(unique(SYMBOL)), collapse = ";"),
      .groups = "drop"
    ) |>
    dplyr::select(gene_set_source, gene_set_id, pathway_name, pathway, linked_genes_n, linked_genes)

  list(pathways = stats::setNames(go_list$genes, go_list$pathway), audit = audit)
}

custom_gene_sets <- list(
  mitophagy_core = c(
    "PINK1", "PRKN", "PARK7", "BNIP3", "BNIP3L", "FUNDC1", "BCL2L13",
    "SQSTM1", "OPTN", "CALCOCO2", "TBK1", "MFN1", "MFN2", "VDAC1",
    "TOMM20", "TOMM22", "TOMM40", "MAP1LC3A", "MAP1LC3B", "GABARAP",
    "GABARAPL1", "GABARAPL2", "ATG5", "ATG7", "ATG12", "ATG13",
    "ATG14", "ULK1", "ULK2", "BECN1", "WIPI1", "WIPI2"
  ),
  autophagy_core = c(
    "ULK1", "ULK2", "ATG13", "RB1CC1", "BECN1", "PIK3C3", "PIK3R4",
    "ATG14", "WIPI1", "WIPI2", "ATG2A", "ATG2B", "ATG5", "ATG7",
    "ATG10", "ATG12", "ATG16L1", "MAP1LC3A", "MAP1LC3B", "GABARAP",
    "GABARAPL1", "GABARAPL2", "SQSTM1", "NBR1", "OPTN", "CALCOCO2"
  ),
  proteasome_core = c(
    paste0("PSMA", 1:7), paste0("PSMB", 1:10), paste0("PSMC", 1:6),
    paste0("PSMD", 1:14), paste0("PSME", 1:4)
  ),
  mitochondrial_oxphos_core = c(
    "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5",
    "MT-ND6", "MT-CYB", "MT-CO1", "MT-CO2", "MT-CO3", "MT-ATP6",
    "MT-ATP8", "NDUFA1", "NDUFA2", "NDUFA3", "NDUFA4", "NDUFA5",
    "NDUFA6", "NDUFA7", "NDUFA8", "NDUFA9", "NDUFA10", "NDUFB1",
    "NDUFB2", "NDUFB3", "NDUFB4", "NDUFB5", "NDUFB6", "NDUFB7",
    "NDUFB8", "NDUFB9", "NDUFB10", "NDUFS1", "NDUFS2", "NDUFS3",
    "NDUFS4", "NDUFS5", "NDUFS6", "NDUFS7", "NDUFS8", "UQCRC1",
    "UQCRC2", "UQCRB", "UQCRQ", "COX4I1", "COX5A", "COX5B",
    "COX6A1", "COX6B1", "COX6C", "COX7A2", "COX7B", "COX7C",
    "COX8A", "ATP5F1A", "ATP5F1B", "ATP5F1C", "ATP5F1D", "ATP5F1E",
    "ATP5MC1", "ATP5MC2", "ATP5MC3", "SDHA", "SDHB", "SDHC", "SDHD"
  ),
  mitochondrial_stress_dynamics = c(
    "MFN1", "MFN2", "OPA1", "DNM1L", "FIS1", "MFF", "MIEF1", "MIEF2",
    "PINK1", "PRKN", "HSPA9", "LONP1", "CLPP", "AFG3L2", "SPG7",
    "ATF4", "ATF5", "DDIT3", "HSPD1", "HSPE1", "HSPA60", "SOD1",
    "SOD2", "GPX1", "GPX4", "PRDX3", "TXN2"
  ),
  mitochondrial_translation = c(
    "TUFM", "TSFM", "GFM1", "GFM2", "MTIF2", "MTIF3", "MTERF1",
    "MRPL1", "MRPL2", "MRPL3", "MRPL4", "MRPL9", "MRPL11", "MRPL12",
    "MRPL13", "MRPL14", "MRPL15", "MRPL16", "MRPL17", "MRPL18",
    "MRPL19", "MRPL20", "MRPL21", "MRPL22", "MRPL23", "MRPL24",
    "MRPL27", "MRPL28", "MRPL30", "MRPL32", "MRPL33", "MRPL34",
    "MRPL35", "MRPL36", "MRPL37", "MRPL38", "MRPL39", "MRPL40",
    "MRPL41", "MRPL42", "MRPL43", "MRPL44", "MRPL45", "MRPL46",
    "MRPL47", "MRPL48", "MRPL49", "MRPL50", "MRPL51", "MRPL52",
    "MRPL53", "MRPL54", "MRPL55", "MRPS2", "MRPS5", "MRPS6",
    "MRPS7", "MRPS9", "MRPS10", "MRPS11", "MRPS12", "MRPS14",
    "MRPS15", "MRPS16", "MRPS17", "MRPS18A", "MRPS18B", "MRPS18C",
    "MRPS21", "MRPS22", "MRPS23", "MRPS24", "MRPS25", "MRPS26",
    "MRPS27", "MRPS28", "MRPS30", "MRPS31", "MRPS33", "MRPS34",
    "MRPS35", "MRPS36"
  )
)
custom_gene_sets <- lapply(custom_gene_sets, unique)

run_fgsea_with_audit <- function(pathway_list_all, audit, rank_stats, source, min_size = 5, max_size = Inf) {
  ranked_genes <- names(rank_stats)
  gene_counts <- vapply(pathway_list_all, function(gs) length(intersect(gs, ranked_genes)), integer(1))
  pathway_sizes <- vapply(pathway_list_all, length, integer(1))
  max_size_for_test <- if (is.infinite(max_size)) max(gene_counts, na.rm = TRUE) else max_size
  testable <- names(gene_counts)[gene_counts >= min_size & gene_counts <= max_size_for_test]

  audit_out <- audit |>
    dplyr::mutate(
      gene_set_source = source,
      ranked_genes_n = gene_counts[pathway],
      ranked_genes_n = ifelse(is.na(ranked_genes_n), 0L, ranked_genes_n),
      gene_set_size_total = pathway_sizes[pathway],
      gene_set_size_total = ifelse(is.na(gene_set_size_total), 0L, gene_set_size_total),
      tested = pathway %in% testable,
      untested_reason = dplyr::case_when(
        tested ~ "",
        ranked_genes_n < min_size ~ paste0("fewer_than_", min_size, "_ranked_genes"),
        ranked_genes_n > max_size_for_test ~ paste0("more_than_", max_size_for_test, "_ranked_genes"),
        TRUE ~ "not_tested"
      )
    )

  if (!length(testable)) {
    return(audit_out |>
             dplyr::mutate(
               size = NA_integer_, ES = NA_real_, NES = NA_real_, pval = NA_real_, padj = NA_real_,
               log2err = NA_real_, direction = "untested", leadingEdge = ""
             ))
  }

  test_list <- lapply(pathway_list_all[testable], function(gs) intersect(gs, ranked_genes))
  fg <- suppressWarnings(fgsea::fgseaMultilevel(
    pathways = test_list,
    stats = rank_stats,
    minSize = min_size,
    maxSize = max_size_for_test
  ))
  fg <- as.data.frame(fg)
  fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  fg <- fg |>
    dplyr::mutate(
      direction = dplyr::case_when(
        is.na(NES) ~ "untested",
        NES > 0 ~ "relapse_up",
        NES < 0 ~ "relapse_down",
        TRUE ~ "flat"
      )
    )

  audit_out |>
    dplyr::left_join(fg, by = "pathway") |>
    dplyr::mutate(
      direction = dplyr::case_when(
        tested & is.na(direction) ~ "tested_no_result",
        !tested ~ "untested",
        TRUE ~ direction
      ),
      leadingEdge = ifelse(is.na(leadingEdge), "", leadingEdge)
    ) |>
    dplyr::arrange(dplyr::desc(tested), padj, dplyr::desc(abs(NES)), pathway)
}

paired_metric_stats <- function(df, metric_col, feature_name, metric_label) {
  wide <- df |>
    dplyr::select(patient_pair, condition_simple, value = dplyr::all_of(metric_col)) |>
    tidyr::pivot_wider(names_from = condition_simple, values_from = value) |>
    dplyr::filter(!is.na(diagnosis), !is.na(relapse))
  if (nrow(wide) < 2) {
    return(tibble::tibble(
      feature = feature_name, metric = metric_label, n_pairs = nrow(wide),
      mean_diagnosis = NA_real_, mean_relapse = NA_real_, mean_delta = NA_real_,
      n_pairs_relapse_higher = NA_integer_, wilcox_p = NA_real_, paired_t_p = NA_real_
    ))
  }
  tibble::tibble(
    feature = feature_name,
    metric = metric_label,
    n_pairs = nrow(wide),
    mean_diagnosis = mean(wide$diagnosis),
    mean_relapse = mean(wide$relapse),
    mean_delta = mean(wide$relapse - wide$diagnosis),
    n_pairs_relapse_higher = sum(wide$relapse > wide$diagnosis),
    wilcox_p = suppressWarnings(stats::wilcox.test(wide$relapse, wide$diagnosis, paired = TRUE, exact = FALSE)$p.value),
    paired_t_p = suppressWarnings(stats::t.test(wide$relapse, wide$diagnosis, paired = TRUE)$p.value)
  )
}

make_rank_for_compartment <- function(compartment_name, cells, counts_all, sample_info, min_cells) {
  meta_sub <- meta_all[match(cells, meta_all$cell), , drop = FALSE]
  count_tbl <- meta_sub |>
    dplyr::count(patient_pair, sample_id, condition_simple, name = "n_cells") |>
    dplyr::right_join(sample_info, by = c("patient_pair", "sample_id", "condition_simple")) |>
    dplyr::mutate(n_cells = ifelse(is.na(n_cells), 0L, n_cells)) |>
    dplyr::arrange(patient_pair, condition_simple)

  valid <- all(count_tbl$n_cells >= min_cells) &&
    length(unique(count_tbl$sample_id[count_tbl$n_cells > 0])) == nrow(sample_info) &&
    all(tapply(count_tbl$n_cells > 0, count_tbl$patient_pair, sum) == 2)
  reason <- if (valid) "" else paste0("requires_all_6_paired_samples_with_at_least_", min_cells, "_cells")

  validity <- tibble::tibble(
    compartment = compartment_name,
    n_cells_total = length(cells),
    n_samples_with_cells = sum(count_tbl$n_cells > 0),
    min_cells_per_sample_observed = min(count_tbl$n_cells),
    median_cells_per_sample_observed = stats::median(count_tbl$n_cells),
    valid_for_paired_pseudobulk = valid,
    skip_reason = reason
  )

  if (!valid) {
    return(list(validity = validity, counts = count_tbl, deg = NULL, rank_stats = NULL))
  }

  sample_cells <- split(meta_sub$cell, meta_sub$sample_id)
  pb_counts <- do.call(cbind, lapply(sample_info$sample_id, function(sid) {
    Matrix::rowSums(counts_all[, sample_cells[[sid]], drop = FALSE])
  }))
  colnames(pb_counts) <- sample_info$sample_id
  rownames(sample_info) <- sample_info$sample_id

  group <- edgeR::DGEList(counts = pb_counts, samples = as.data.frame(sample_info))
  design <- stats::model.matrix(~ patient_pair + condition_simple, data = sample_info)
  keep <- edgeR::filterByExpr(group, design = design)
  if (sum(keep) < 100) {
    validity$valid_for_paired_pseudobulk <- FALSE
    validity$skip_reason <- "fewer_than_100_genes_after_filterByExpr"
    return(list(validity = validity, counts = count_tbl, deg = NULL, rank_stats = NULL))
  }
  group <- group[keep, , keep.lib.sizes = FALSE]
  group <- edgeR::calcNormFactors(group)
  group <- edgeR::estimateDisp(group, design, robust = TRUE)
  fit <- edgeR::glmQLFit(group, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = "condition_simplerelapse")
  deg <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
  deg$gene <- rownames(deg)
  deg$FDR <- stats::p.adjust(deg$PValue, method = "BH")
  deg$rank_stat <- sign(deg$logFC) * sqrt(pmax(deg$F, 0))
  deg <- deg |>
    dplyr::select(gene, logFC, logCPM, F, PValue, FDR, rank_stat) |>
    dplyr::arrange(dplyr::desc(rank_stat))
  rank_df <- deg |>
    dplyr::filter(is.finite(rank_stat), !duplicated(gene))
  rank_stats <- rank_df$rank_stat
  names(rank_stats) <- rank_df$gene
  rank_stats <- sort(rank_stats, decreasing = TRUE)

  list(validity = validity, counts = count_tbl, deg = deg, rank_stats = rank_stats)
}

compute_custom_module_stats <- function(compartment_name, cells, data_all, counts_all, sample_info, gene_sets, cell_total_counts) {
  meta_sub <- meta_all[match(cells, meta_all$cell), , drop = FALSE]
  sample_cells <- split(meta_sub$cell, meta_sub$sample_id)
  rows <- list()
  for (set_name in names(gene_sets)) {
    genes <- intersect(gene_sets[[set_name]], rownames(data_all))
    if (!length(genes)) next
    sample_rows <- lapply(sample_info$sample_id, function(sid) {
      if (!sid %in% names(sample_cells) || !length(sample_cells[[sid]])) {
        return(tibble::tibble(
          compartment = compartment_name, gene_set = set_name, sample_id = sid,
          mean_lognorm_module = NA_real_, pct_cells_any_detected = NA_real_,
          pseudobulk_mean_logCPM = NA_real_, genes_present_n = length(genes),
          genes_present = paste(genes, collapse = ";")
        ))
      }
      sub_cells <- sample_cells[[sid]]
      cell_scores <- Matrix::colMeans(data_all[genes, sub_cells, drop = FALSE])
      detected_any <- Matrix::colSums(counts_all[genes, sub_cells, drop = FALSE] > 0) > 0
      pb_gene_counts <- Matrix::rowSums(counts_all[genes, sub_cells, drop = FALSE])
      lib_size <- sum(cell_total_counts[sub_cells])
      tibble::tibble(
        compartment = compartment_name,
        gene_set = set_name,
        sample_id = sid,
        mean_lognorm_module = mean(cell_scores, na.rm = TRUE),
        pct_cells_any_detected = mean(detected_any, na.rm = TRUE) * 100,
        pseudobulk_mean_logCPM = mean(log2((as.numeric(pb_gene_counts) + 0.5) / (lib_size + 1) * 1e6), na.rm = TRUE),
        genes_present_n = length(genes),
        genes_present = paste(genes, collapse = ";")
      )
    })
    sample_df <- dplyr::bind_rows(sample_rows) |>
      dplyr::left_join(sample_info, by = "sample_id")
    stats_df <- dplyr::bind_rows(
      paired_metric_stats(sample_df, "mean_lognorm_module", set_name, "mean_lognorm_module"),
      paired_metric_stats(sample_df, "pct_cells_any_detected", set_name, "pct_cells_any_detected"),
      paired_metric_stats(sample_df, "pseudobulk_mean_logCPM", set_name, "pseudobulk_mean_logCPM")
    ) |>
      dplyr::mutate(compartment = compartment_name, .before = 1)
    rows[[set_name]] <- list(sample = sample_df, stats = stats_df)
  }
  list(
    sample = dplyr::bind_rows(lapply(rows, `[[`, "sample")),
    stats = dplyr::bind_rows(lapply(rows, `[[`, "stats"))
  )
}

log_message("Reading annotated Seurat object")
obj_file <- file.path(paths$objects, "GSE262271_seurat_qc_harmony_annotated.rds")
if (!file.exists(obj_file)) stop("Missing annotated Seurat object: ", obj_file, call. = FALSE)
obj <- readRDS(obj_file)
obj <- join_layers_safe(obj)
DefaultAssay(obj) <- "RNA"
if (!"cluster_res0.6" %in% colnames(obj@meta.data) && "RNA_snn_res.0.6" %in% colnames(obj@meta.data)) {
  obj$cluster_res0.6 <- as.character(obj$RNA_snn_res.0.6)
}

meta_all <- tibble::rownames_to_column(obj@meta.data, "cell")
required_meta <- c("sample_id", "patient_pair", "condition_simple", "cluster_res0.6", "major_annotation", "broad_annotation")
missing_meta <- setdiff(required_meta, colnames(meta_all))
if (length(missing_meta)) stop("Annotated object is missing metadata columns: ", paste(missing_meta, collapse = ", "), call. = FALSE)
meta_all$condition_simple <- as.character(meta_all$condition_simple)
meta_all$condition_simple <- factor(meta_all$condition_simple, levels = c("diagnosis", "relapse"))
meta_all$patient_pair <- factor(meta_all$patient_pair)
meta_all$cluster_res0.6 <- as.character(meta_all$cluster_res0.6)

sample_info <- meta_all |>
  dplyr::distinct(sample_id, patient_pair, condition_simple) |>
  dplyr::filter(!is.na(condition_simple)) |>
  dplyr::arrange(patient_pair, condition_simple)
if (nrow(sample_info) != 6) {
  stop("Expected 6 paired samples for GSE262271; found ", nrow(sample_info), ".", call. = FALSE)
}
if (any(duplicated(sample_info$sample_id))) stop("Duplicated sample_id in sample metadata.", call. = FALSE)

counts_all <- get_assay_matrix(obj, "RNA", "counts")
data_all <- get_assay_matrix(obj, "RNA", "data")
cell_total_counts <- Matrix::colSums(counts_all)

cnv_file <- file.path(paths$results, "infercnv_nonTref_cell_cnv_burden.csv")
if (file.exists(cnv_file)) {
  cnv <- utils::read.csv(cnv_file, stringsAsFactors = FALSE, check.names = FALSE)
  cell_col <- intersect(c("cell", "barcode", "cell_id", "cell_name"), colnames(cnv))
  if (!length(cell_col)) cell_col <- colnames(cnv)[[1]]
  score_col <- intersect(c("infercnv_mean_abs_log2", "mean_abs_log2"), colnames(cnv))
  if (!length(score_col)) score_col <- grep("mean.*abs.*log2|abs.*log2", colnames(cnv), ignore.case = TRUE, value = TRUE)
  if (!length(score_col)) {
    log_message("inferCNV burden file found but no usable score column detected; CNV-high compartments skipped")
    meta_all$infercnv_mean_abs_log2 <- NA_real_
  } else {
    cnv_score <- cnv[[score_col[[1]]]]
    names(cnv_score) <- cnv[[cell_col[[1]]]]
    meta_all$infercnv_mean_abs_log2 <- as.numeric(cnv_score[meta_all$cell])
  }
} else {
  log_message("inferCNV burden file not found; CNV-high compartments skipped:", cnv_file)
  meta_all$infercnv_mean_abs_log2 <- NA_real_
}

candidate <- meta_all$broad_annotation == "T_lineage_candidate_malignant"
compartment_cells <- list(
  broad_T_lineage_candidate_malignant = meta_all$cell[candidate]
)
candidate_cnv <- meta_all$infercnv_mean_abs_log2[candidate]
if (sum(is.finite(candidate_cnv)) > 100) {
  q50 <- stats::quantile(candidate_cnv, 0.50, na.rm = TRUE)
  q75 <- stats::quantile(candidate_cnv, 0.75, na.rm = TRUE)
  compartment_cells$CNV_high_top50_candidate <- meta_all$cell[candidate & is.finite(meta_all$infercnv_mean_abs_log2) & meta_all$infercnv_mean_abs_log2 >= q50]
  compartment_cells$CNV_high_top25_candidate <- meta_all$cell[candidate & is.finite(meta_all$infercnv_mean_abs_log2) & meta_all$infercnv_mean_abs_log2 >= q75]
}

major_targets <- c(
  "Immature_HSC_like_TALL",
  "Cycling_TALL_blast_like",
  "Quiescent_CD44_high_TALL_like",
  "ER_stress_TALL_like",
  "TALL_blast_like",
  "T_lineage_like"
)
for (target in intersect(major_targets, unique(meta_all$major_annotation))) {
  compartment_cells[[target]] <- meta_all$cell[meta_all$major_annotation == target]
}

t_clusters <- sort(unique(meta_all$cluster_res0.6[candidate]))
for (cl in t_clusters) {
  cl_major <- meta_all |>
    dplyr::filter(cluster_res0.6 == cl) |>
    dplyr::count(major_annotation, sort = TRUE) |>
    dplyr::slice_head(n = 1) |>
    dplyr::pull(major_annotation)
  nm <- paste0("cluster_res0.6_", cl, "_", slugify(cl_major))
  compartment_cells[[nm]] <- meta_all$cell[meta_all$cluster_res0.6 == cl & candidate]
}
compartment_cells <- compartment_cells[lengths(compartment_cells) > 0]

log_message("Building gene sets")
kegg_sets <- build_kegg_gene_sets(paths)
go_bp_sets <- build_go_bp_gene_sets(rownames(counts_all))
custom_audit <- tibble::tibble(
  gene_set_source = "custom_mito",
  gene_set_id = names(custom_gene_sets),
  pathway_name = names(custom_gene_sets),
  pathway = names(custom_gene_sets),
  linked_genes_n = lengths(custom_gene_sets),
  linked_genes = vapply(custom_gene_sets, paste, collapse = ";", FUN.VALUE = character(1))
)

focus_kegg_pattern <- paste(c(
  "mitophagy", "autophagy", "oxidative phosphorylation", "proteasome",
  "ubiquitin", "thermogenesis", "reactive oxygen", "glutathione", "peroxisome",
  "tca cycle", "citrate cycle", "pyruvate metabolism", "jak-stat", "hif-1"
), collapse = "|")
focus_go_pattern <- paste(c(
  "mitophagy", "mitochond", "autophagy", "respiratory chain", "oxidative phosphorylation",
  "electron transport", "proteasome", "ubiquitin", "reactive oxygen", "cellular respiration",
  "atp synthesis", "tricarboxylic", "pyruvate", "apoptotic", "jak-stat", "hypoxia"
), collapse = "|")

all_counts <- list()
all_validity <- list()
all_deg <- list()
all_kegg <- list()
all_go <- list()
all_custom_gsea <- list()
all_custom_sample <- list()
all_custom_stats <- list()

for (compartment_name in names(compartment_cells)) {
  log_message("Running sensitivity compartment:", compartment_name)
  cells <- compartment_cells[[compartment_name]]
  comp <- make_rank_for_compartment(compartment_name, cells, counts_all, sample_info, min_cells_per_sample)
  all_counts[[compartment_name]] <- comp$counts |>
    dplyr::mutate(compartment = compartment_name, .before = 1)
  all_validity[[compartment_name]] <- comp$validity

  module_out <- compute_custom_module_stats(compartment_name, cells, data_all, counts_all, sample_info, custom_gene_sets, cell_total_counts)
  all_custom_sample[[compartment_name]] <- module_out$sample
  all_custom_stats[[compartment_name]] <- module_out$stats

  if (is.null(comp$rank_stats)) next

  deg_file <- file.path(paths$sensitivity_results, paste0("DEG_pseudobulk_", slugify(compartment_name), "_relapse_vs_diagnosis.csv"))
  write_csv_file(dplyr::mutate(comp$deg, compartment = compartment_name, .before = 1), deg_file)
  all_deg[[compartment_name]] <- dplyr::mutate(comp$deg, compartment = compartment_name, .before = 1)

  kegg_res <- run_fgsea_with_audit(kegg_sets$pathways, kegg_sets$audit, comp$rank_stats, "KEGG", min_size = 5, max_size = Inf) |>
    dplyr::mutate(compartment = compartment_name, .before = 1)
  all_kegg[[compartment_name]] <- kegg_res

  go_res <- run_fgsea_with_audit(go_bp_sets$pathways, go_bp_sets$audit, comp$rank_stats, "GO_BP", min_size = 10, max_size = 500) |>
    dplyr::mutate(compartment = compartment_name, .before = 1)
  all_go[[compartment_name]] <- go_res

  custom_present <- lapply(custom_gene_sets, function(gs) intersect(gs, names(comp$rank_stats)))
  custom_present <- custom_present[lengths(custom_present) >= 3]
  if (length(custom_present)) {
    custom_res <- run_fgsea_with_audit(custom_present, custom_audit, comp$rank_stats, "custom_mito", min_size = 3, max_size = Inf) |>
      dplyr::mutate(compartment = compartment_name, .before = 1)
    all_custom_gsea[[compartment_name]] <- custom_res
  }
}

cell_counts <- dplyr::bind_rows(all_counts)
validity <- dplyr::bind_rows(all_validity)
deg_all <- dplyr::bind_rows(all_deg)
kegg_all <- dplyr::bind_rows(all_kegg)
go_all <- dplyr::bind_rows(all_go)
custom_gsea_all <- dplyr::bind_rows(all_custom_gsea)
custom_sample_all <- dplyr::bind_rows(all_custom_sample)
custom_stats_all <- dplyr::bind_rows(all_custom_stats)

write_csv_file(cell_counts, file.path(paths$sensitivity_results, "sensitivity_compartment_cell_counts.csv"))
write_csv_file(validity, file.path(paths$sensitivity_results, "sensitivity_compartment_validity.csv"))
write_csv_file(deg_all, file.path(paths$sensitivity_results, "sensitivity_all_compartments_pseudobulk_DEG_rankings.csv"))
write_csv_file(kegg_all, file.path(paths$sensitivity_results, "sensitivity_kegg_all_pathways_by_compartment.csv"))
write_csv_file(
  kegg_all |>
    dplyr::filter(grepl(focus_kegg_pattern, pathway_name, ignore.case = TRUE)) |>
    dplyr::arrange(compartment, dplyr::desc(tested), padj, dplyr::desc(abs(NES))),
  file.path(paths$sensitivity_results, "sensitivity_kegg_focus_pathways_by_compartment.csv")
)
write_csv_file(
  go_all |>
    dplyr::select(-dplyr::any_of("linked_genes")),
  file.path(paths$sensitivity_results, "sensitivity_go_bp_all_terms_by_compartment.csv")
)
write_csv_file(
  go_all |>
    dplyr::filter(grepl(focus_go_pattern, pathway_name, ignore.case = TRUE)) |>
    dplyr::arrange(compartment, dplyr::desc(tested), padj, dplyr::desc(abs(NES))),
  file.path(paths$sensitivity_results, "sensitivity_go_bp_focus_terms_by_compartment.csv")
)
write_csv_file(custom_gsea_all, file.path(paths$sensitivity_results, "sensitivity_custom_mito_gene_sets_fgsea_by_compartment.csv"))
write_csv_file(custom_sample_all, file.path(paths$sensitivity_results, "sensitivity_custom_mito_module_sample_scores_by_compartment.csv"))
write_csv_file(custom_stats_all, file.path(paths$sensitivity_results, "sensitivity_custom_mito_module_paired_statistics_by_compartment.csv"))

summary_focus <- dplyr::bind_rows(
  kegg_all |>
    dplyr::filter(grepl("mitophagy|oxidative phosphorylation|proteasome|autophagy|reactive oxygen|glutathione|thermogenesis", pathway_name, ignore.case = TRUE)) |>
    dplyr::transmute(
      source = "KEGG", compartment, gene_set = pathway_name, tested, NES, padj, direction,
      leadingEdge, ranked_genes_n, note = ""
    ),
  go_all |>
    dplyr::filter(grepl("mitophagy|mitochond|autophagy|respiratory chain|oxidative phosphorylation|electron transport|proteasome|ubiquitin|reactive oxygen", pathway_name, ignore.case = TRUE)) |>
    dplyr::transmute(
      source = "GO_BP", compartment, gene_set = pathway_name, tested, NES, padj, direction,
      leadingEdge, ranked_genes_n, note = ""
    ),
  custom_gsea_all |>
    dplyr::transmute(
      source = "custom", compartment, gene_set = pathway_name, tested, NES, padj, direction,
      leadingEdge, ranked_genes_n, note = ""
    )
) |>
  dplyr::arrange(source, gene_set, compartment)
write_csv_file(summary_focus, file.path(paths$sensitivity_results, "sensitivity_summary_mitophagy_mito_proteasome.csv"))

valid_compartments <- validity |>
  dplyr::filter(valid_for_paired_pseudobulk) |>
  dplyr::pull(compartment)

if (nrow(kegg_all)) {
  kegg_plot <- kegg_all |>
    dplyr::filter(
      compartment %in% valid_compartments,
      grepl("Mitophagy|Oxidative phosphorylation|Proteasome|Autophagy|Reactive oxygen|Glutathione|Thermogenesis", pathway_name, ignore.case = TRUE),
      tested
    ) |>
    dplyr::mutate(
      pathway_name = factor(pathway_name, levels = rev(sort(unique(pathway_name)))),
      compartment = factor(compartment, levels = valid_compartments),
      neg_log10_fdr = -log10(pmax(padj, 1e-300))
    )
  if (nrow(kegg_plot)) {
    p <- ggplot2::ggplot(kegg_plot, ggplot2::aes(x = compartment, y = pathway_name, fill = NES, size = neg_log10_fdr)) +
      ggplot2::geom_point(shape = 21, color = "grey25", stroke = 0.2) +
      ggplot2::scale_fill_gradient2(low = "#356A8A", mid = "white", high = "#B2433F", midpoint = 0) +
      ggplot2::scale_size_continuous(range = c(1.2, 5.0)) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = NULL, fill = "NES", size = "-log10 FDR") +
      theme_pub(8) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    ggplot2::ggsave(file.path(paths$sensitivity_figures, "sensitivity_kegg_focus_dotplot.png"), p, width = 10, height = max(4.5, 0.28 * length(unique(kegg_plot$pathway_name))), dpi = 300)
  }
}

if (nrow(custom_gsea_all)) {
  custom_plot <- custom_gsea_all |>
    dplyr::filter(compartment %in% valid_compartments, tested) |>
    dplyr::mutate(
      pathway_name = factor(pathway_name, levels = rev(names(custom_gene_sets))),
      compartment = factor(compartment, levels = valid_compartments),
      neg_log10_fdr = -log10(pmax(padj, 1e-300))
    )
  if (nrow(custom_plot)) {
    p <- ggplot2::ggplot(custom_plot, ggplot2::aes(x = compartment, y = pathway_name, fill = NES, size = neg_log10_fdr)) +
      ggplot2::geom_point(shape = 21, color = "grey25", stroke = 0.2) +
      ggplot2::scale_fill_gradient2(low = "#356A8A", mid = "white", high = "#B2433F", midpoint = 0) +
      ggplot2::scale_size_continuous(range = c(1.2, 5.0)) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = NULL, fill = "NES", size = "-log10 FDR") +
      theme_pub(8) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    ggplot2::ggsave(file.path(paths$sensitivity_figures, "sensitivity_custom_mito_gene_sets_dotplot.png"), p, width = 9.5, height = 4.8, dpi = 300)
  }
}

if (nrow(go_all)) {
  go_focus_plot <- go_all |>
    dplyr::filter(
      compartment %in% valid_compartments,
      grepl("mitophagy|mitochondrial respiratory chain|oxidative phosphorylation|electron transport chain|proteasome|macroautophagy", pathway_name, ignore.case = TRUE),
      tested
    ) |>
    dplyr::group_by(compartment) |>
    dplyr::arrange(padj, dplyr::desc(abs(NES)), .by_group = TRUE) |>
    dplyr::slice_head(n = 15) |>
    dplyr::ungroup() |>
    dplyr::distinct(compartment, pathway_name, .keep_all = TRUE) |>
    dplyr::mutate(
      pathway_name = factor(pathway_name, levels = rev(unique(pathway_name[order(NES)]))),
      compartment = factor(compartment, levels = valid_compartments),
      neg_log10_fdr = -log10(pmax(padj, 1e-300))
    )
  if (nrow(go_focus_plot)) {
    p <- ggplot2::ggplot(go_focus_plot, ggplot2::aes(x = compartment, y = pathway_name, fill = NES, size = neg_log10_fdr)) +
      ggplot2::geom_point(shape = 21, color = "grey25", stroke = 0.2) +
      ggplot2::scale_fill_gradient2(low = "#356A8A", mid = "white", high = "#B2433F", midpoint = 0) +
      ggplot2::scale_size_continuous(range = c(1.0, 4.5)) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = NULL, fill = "NES", size = "-log10 FDR") +
      theme_pub(7) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    ggplot2::ggsave(file.path(paths$sensitivity_figures, "sensitivity_go_bp_mito_focus_dotplot.png"), p, width = 11, height = max(5.5, 0.18 * nrow(go_focus_plot)), dpi = 300)
  }
}

sig_kegg <- kegg_all |>
  dplyr::filter(tested, !is.na(padj), padj < 0.05) |>
  dplyr::arrange(compartment, padj)
sig_go <- go_all |>
  dplyr::filter(tested, !is.na(padj), padj < 0.05) |>
  dplyr::arrange(compartment, padj)
mitophagy_rows <- summary_focus |>
  dplyr::filter(grepl("mitophagy|mitophagy_core", gene_set, ignore.case = TRUE)) |>
  dplyr::arrange(source, compartment)

report_lines <- c(
  "# Sensitivity analysis: relapse vs diagnosis in candidate malignant-like T-ALL compartments",
  "",
  "## Design",
  paste0("- Dataset: GSE262271, paired diagnosis/relapse T-ALL samples, N = ", length(unique(sample_info$patient_pair)), " patient pairs."),
  paste0("- Minimum cells per sample required for paired pseudobulk DEG/GSEA: ", min_cells_per_sample, "."),
  "- Annotation sensitivity compartments include broad T-lineage candidate malignant cells, inferCNV-high candidate cells when available, major marker-panel states, and individual res0.6 T candidate clusters.",
  "- Differential expression ranking: edgeR paired pseudobulk quasi-likelihood model `~ patient_pair + condition_simple`; rank statistic is `sign(logFC) * sqrt(F)`.",
  "- Gene set tests: fgsea over KEGG, GO biological process, and curated custom mitophagy/mitochondrial/proteasome gene sets.",
  "",
  "## Compartment validity",
  paste0("- Compartments evaluated: ", nrow(validity), "."),
  paste0("- Compartments passing paired pseudobulk criteria: ", sum(validity$valid_for_paired_pseudobulk), "."),
  paste0("- Valid compartments: ", paste(valid_compartments, collapse = "; ")),
  "",
  "## Significant pathway counts",
  paste0("- KEGG significant rows across valid compartments (FDR < 0.05): ", nrow(sig_kegg), "."),
  paste0("- GO BP significant rows across valid compartments (FDR < 0.05): ", nrow(sig_go), "."),
  "",
  "## Mitophagy-specific rows",
  if (nrow(mitophagy_rows)) {
    paste0(
      "- ", mitophagy_rows$source, " / ", mitophagy_rows$compartment, " / ", mitophagy_rows$gene_set,
      ": NES ", signif(mitophagy_rows$NES, 4),
      ", FDR ", signif(mitophagy_rows$padj, 4),
      ", ", mitophagy_rows$direction, "."
    )
  } else {
    "- No mitophagy-focused row was available after gene-set size filtering."
  },
  "",
  "## Output files",
  "- results/sensitivity_mitophagy_gsea/sensitivity_compartment_cell_counts.csv",
  "- results/sensitivity_mitophagy_gsea/sensitivity_compartment_validity.csv",
  "- results/sensitivity_mitophagy_gsea/sensitivity_all_compartments_pseudobulk_DEG_rankings.csv",
  "- results/sensitivity_mitophagy_gsea/sensitivity_kegg_all_pathways_by_compartment.csv",
  "- results/sensitivity_mitophagy_gsea/sensitivity_kegg_focus_pathways_by_compartment.csv",
  "- results/sensitivity_mitophagy_gsea/sensitivity_go_bp_all_terms_by_compartment.csv",
  "- results/sensitivity_mitophagy_gsea/sensitivity_go_bp_focus_terms_by_compartment.csv",
  "- results/sensitivity_mitophagy_gsea/sensitivity_custom_mito_gene_sets_fgsea_by_compartment.csv",
  "- results/sensitivity_mitophagy_gsea/sensitivity_custom_mito_module_sample_scores_by_compartment.csv",
  "- results/sensitivity_mitophagy_gsea/sensitivity_custom_mito_module_paired_statistics_by_compartment.csv",
  "- results/sensitivity_mitophagy_gsea/sensitivity_summary_mitophagy_mito_proteasome.csv",
  "- figures/sensitivity_mitophagy_gsea/sensitivity_kegg_focus_dotplot.png",
  "- figures/sensitivity_mitophagy_gsea/sensitivity_go_bp_mito_focus_dotplot.png",
  "- figures/sensitivity_mitophagy_gsea/sensitivity_custom_mito_gene_sets_dotplot.png"
)
writeLines(report_lines, file.path(paths$sensitivity_results, "SENSITIVITY_MITOPHAGY_GSEA_REPORT.md"), useBytes = TRUE)

capture.output(sessionInfo(), file = file.path(paths$logs, paste0("sessionInfo_sensitivity_mitophagy_gsea_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")))
log_message("Sensitivity mitophagy/GSEA analysis complete")
