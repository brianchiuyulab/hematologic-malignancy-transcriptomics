options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args) >= 1) args[[1]] else "outputs/GSE262271_TALL_relapse_scRNAseq_final", winslash = "/", mustWork = TRUE)
bulk_root <- normalizePath(if (length(args) >= 2) args[[2]] else Sys.getenv("GSE218858_BULK_ROOT", unset = "GSE218858_R_DESeq2_final_deliverables"), winslash = "/", mustWork = TRUE)
lib <- Sys.getenv("R_LIBS_USER", unset = "C:/Users/User/Documents/R/win-library/4.4")
if (dir.exists(lib)) .libPaths(unique(c(normalizePath(lib, winslash = "/"), .libPaths())))
for (p in c("dplyr", "ggplot2", "patchwork")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)

csv_dir <- file.path(root, "csv", "potential_malignant_T_cells")
png_dir <- file.path(root, "png", "potential_malignant_T_cells")
sc <- read.csv(file.path(csv_dir, "KEGG_GSEA_all_pathways_relapse_vs_diagnosis.csv"), check.names = FALSE)
bulk <- read.csv(file.path(bulk_root, "csv_results", "gsea", "R_fgsea_all_pathways_Hallmark_KEGG_GO_BP.csv"), check.names = FALSE)
bulk <- bulk[bulk$collection == "KEGG_2019_Mouse_Enrichr_mapped", , drop = FALSE]

canonical_name <- function(x) {
  x <- trimws(x)
  x[x == "Mitophagy - animal"] <- "Mitophagy"
  x[x == "Autophagy - animal"] <- "Autophagy"
  x
}
key_name <- function(x) gsub("[^a-z0-9]+", "", tolower(canonical_name(x)))

sc2 <- sc |>
  dplyr::filter(tested %in% c(TRUE, "TRUE"), is.finite(padj), padj < 0.10) |>
  dplyr::transmute(
    original_sc = pathway_name,
    canonical = canonical_name(pathway_name),
    key = key_name(pathway_name),
    sc_NES = NES,
    sc_FDR = padj,
    sc_direction = ifelse(NES > 0, "up", "down")
  ) |>
  dplyr::group_by(key) |>
  dplyr::slice_min(sc_FDR, n = 1, with_ties = FALSE) |>
  dplyr::ungroup()
bulk2 <- bulk |>
  dplyr::filter(is.finite(padj), padj < 0.10) |>
  dplyr::transmute(
    original_bulk = pathway,
    canonical = canonical_name(pathway),
    key = key_name(pathway),
    bulk_NES = NES,
    bulk_FDR = padj,
    bulk_direction = ifelse(NES > 0, "up", "down")
  ) |>
  dplyr::group_by(key) |>
  dplyr::slice_min(bulk_FDR, n = 1, with_ties = FALSE) |>
  dplyr::ungroup()

membership <- dplyr::full_join(sc2, bulk2, by = "key", suffix = c("_sc", "_bulk")) |>
  dplyr::mutate(
    pathway = dplyr::coalesce(canonical_sc, canonical_bulk),
    in_sc_FDR_lt_0.10 = !is.na(sc_FDR),
    in_bulk_FDR_lt_0.10 = !is.na(bulk_FDR),
    same_direction = in_sc_FDR_lt_0.10 & in_bulk_FDR_lt_0.10 & sc_direction == bulk_direction,
    shared_direction = ifelse(same_direction, sc_direction, NA_character_)
  ) |>
  dplyr::select(pathway, original_sc, original_bulk, in_sc_FDR_lt_0.10, in_bulk_FDR_lt_0.10,
                sc_NES, sc_FDR, sc_direction, bulk_NES, bulk_FDR, bulk_direction,
                same_direction, shared_direction) |>
  dplyr::arrange(dplyr::desc(same_direction), shared_direction, pathway)

write.csv(
  membership,
  file.path(csv_dir, "cross_dataset_KEGG_FDR0.10_pathway_membership.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
shared <- membership |>
  dplyr::filter(same_direction) |>
  dplyr::arrange(shared_direction, pmax(sc_FDR, bulk_FDR), pathway)
write.csv(
  shared,
  file.path(csv_dir, "cross_dataset_KEGG_FDR0.10_same_direction_intersection.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

counts_for <- function(direction) {
  sc_set <- unique(sc2$key[sc2$sc_direction == direction])
  bulk_set <- unique(bulk2$key[bulk2$bulk_direction == direction])
  data.frame(
    direction = direction,
    sc_total = length(sc_set),
    bulk_total = length(bulk_set),
    shared = length(intersect(sc_set, bulk_set)),
    sc_only = length(setdiff(sc_set, bulk_set)),
    bulk_only = length(setdiff(bulk_set, sc_set)),
    stringsAsFactors = FALSE
  )
}
counts <- rbind(counts_for("up"), counts_for("down"))
write.csv(
  counts,
  file.path(csv_dir, "cross_dataset_KEGG_FDR0.10_same_direction_venn_counts.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

circle_data <- function(cx, cy, r, group) {
  theta <- seq(0, 2 * pi, length.out = 361)
  data.frame(x = cx + r * cos(theta), y = cy + r * sin(theta), group = group)
}
venn_panel <- function(z, title) {
  has_overlap <- z$shared > 0
  cx <- if (has_overlap) c(-0.72, 0.72) else c(-1.35, 1.35)
  r <- 1.25
  circles <- rbind(circle_data(cx[[1]], 0, r, "Relapse scRNA-seq"), circle_data(cx[[2]], 0, r, "STAT5B N642H bulk"))
  labels <- data.frame(
    x = c(cx[[1]] - 0.48, 0, cx[[2]] + 0.48), y = c(0, 0, 0),
    label = c(z$sc_only, z$shared, z$bulk_only),
    stringsAsFactors = FALSE
  )
  ggplot2::ggplot(circles, ggplot2::aes(x, y, group = group, fill = group, color = group)) +
    ggplot2::geom_polygon(alpha = 0.20, linewidth = 0.8) +
    ggplot2::geom_text(data = labels, ggplot2::aes(x, y, label = label), inherit.aes = FALSE, fontface = "bold", size = 5) +
    ggplot2::annotate("text", x = cx[[1]] - 0.50, y = 1.48, label = paste0("Relapse scRNA-seq\n(n = ", z$sc_total, ")"), size = 3.5, fontface = "bold") +
    ggplot2::annotate("text", x = cx[[2]] + 0.50, y = 1.48, label = paste0("STAT5B N642H bulk\n(n = ", z$bulk_total, ")"), size = 3.5, fontface = "bold") +
    ggplot2::scale_fill_manual(values = c("Relapse scRNA-seq" = "#3C6E9E", "STAT5B N642H bulk" = "#D07B52")) +
    ggplot2::scale_color_manual(values = c("Relapse scRNA-seq" = "#3C6E9E", "STAT5B N642H bulk" = "#D07B52")) +
    ggplot2::coord_equal(xlim = c(-2.8, 2.8), ylim = c(-1.45, 1.85), clip = "off") +
    ggplot2::labs(title = title) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5), legend.position = "none")
}

p_up <- venn_panel(counts[counts$direction == "up", ], "Enriched in the same upward direction")
p_down <- venn_panel(counts[counts$direction == "down", ], "Enriched in the same downward direction")
combined <- (p_up | p_down) +
  patchwork::plot_annotation(
    title = "Shared KEGG GSEA pathways at FDR < 0.10",
    subtitle = "Relapse potential malignant T cells vs STAT5B N642H bulk RNA-seq",
    caption = "Each dataset is filtered independently at FDR < 0.10; the center contains pathways enriched in the same direction.",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 10),
      plot.caption = ggplot2::element_text(size = 8, hjust = 0)
    )
  )
ggplot2::ggsave(
  file.path(png_dir, "cross_dataset_KEGG_FDR0.10_same_direction_venn.png"),
  combined, width = 10.5, height = 5.5, dpi = 400, bg = "white"
)

print(counts)
cat("\nShared pathways:\n")
print(shared[, c("pathway", "shared_direction", "sc_NES", "sc_FDR", "bulk_NES", "bulk_FDR")], row.names = FALSE)
