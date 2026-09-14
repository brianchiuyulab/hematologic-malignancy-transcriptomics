#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else "C:/Users/User/Documents/BTZ resistance/GSE262271_project_link"
output_root <- if (length(args) >= 2) args[[2]] else file.path(project_root, "_publication_figure_build")

table_root <- file.path(project_root, "02_Tables")
old_figure_root <- file.path(project_root, "01_Figures")
object_file <- file.path(project_root, "04_R_Objects", "GSE262271_seurat_qc_harmony_annotated.rds")

main_dir <- file.path(output_root, "01_Main_Figures")
supp_dir <- file.path(output_root, "02_Supplementary_Figures")
sens_dir <- file.path(output_root, "03_Sensitivity_Checks")
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sens_dir, recursive = TRUE, showWarnings = FALSE)

COL_DX <- "#0072B2"
COL_REL <- "#D55E00"
COL_BULK <- "#E69F00"
COL_UP <- "#D55E00"
COL_DOWN <- "#0072B2"
COL_NEUTRAL <- "#747474"
COL_LIGHT <- "#D9D9D9"

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.25), margin = margin(b = 5)),
      plot.subtitle = element_text(color = "#3F3F3F", size = rel(0.92), margin = margin(b = 10)),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "#202020"),
      legend.title = element_text(face = "bold"),
      legend.key.height = unit(0.48, "cm"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = "#202020"),
      plot.margin = margin(12, 16, 12, 12)
    )
}

save_plot <- function(plot, filename, width, height) {
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 360,
    bg = "white",
    limitsize = FALSE
  )
}

read_table <- function(...) {
  read.csv(file.path(table_root, ...), check.names = FALSE, na.strings = c("", "NA"))
}

fdr_status <- function(q) {
  case_when(
    !is.na(q) & q < 0.05 ~ "FDR < 0.05",
    !is.na(q) & q < 0.10 ~ "Exploratory FDR < 0.10",
    TRUE ~ "Not significant"
  )
}

fdr_mark <- function(q) {
  case_when(
    !is.na(q) & q < 0.05 ~ "*",
    !is.na(q) & q < 0.10 ~ "#",
    TRUE ~ "ns"
  )
}

fmt_q <- function(q) {
  ifelse(is.na(q), "NA", formatC(q, format = "g", digits = 3))
}

status_shapes <- c("FDR < 0.05" = 16, "Exploratory FDR < 0.10" = 17, "Not significant" = 1)

# -----------------------------------------------------------------------------
# Figure 1. Cohort cell counts and paired proportions
# -----------------------------------------------------------------------------
counts <- read_table("01_Key_Results", "potential_malignant_T_cell_counts_and_proportions_by_sample.csv") %>%
  mutate(
    condition = factor(condition_simple, levels = c("diagnosis", "relapse"), labels = c("Diagnosis", "Relapse")),
    sample_label = paste0(patient_pair, "  ", ifelse(condition == "Diagnosis", "Dx", "Rel")),
    other_candidate_T = candidate_T_cells - potential_malignant_T_cells,
    sample_label = factor(sample_label, levels = rev(c("P1  Dx", "P1  Rel", "P2  Dx", "P2  Rel", "P3  Dx", "P3  Rel")))
  )

p_count <- ggplot(counts, aes(y = sample_label)) +
  geom_col(aes(x = candidate_T_cells), width = 0.62, fill = "#E3E3E3", color = NA) +
  geom_col(aes(x = potential_malignant_T_cells, fill = condition), width = 0.62, color = "white", linewidth = 0.25) +
  geom_text(
    aes(x = candidate_T_cells + max(candidate_T_cells) * 0.035,
        label = paste0(comma(potential_malignant_T_cells), " / ", comma(candidate_T_cells), "  (", sprintf("%.1f", potential_malignant_T_cell_percent), "%)")),
    hjust = 0, size = 3.25, family = "Arial"
  ) +
  scale_fill_manual(values = c("Diagnosis" = COL_DX, "Relapse" = COL_REL)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.38)), labels = comma) +
  labs(
    title = "A  Cell counts within the candidate T-cell compartment",
    subtitle = "Colored bars show potential malignant T cells; gray bars show all candidate T-lineage cells",
    x = "Number of cells", y = NULL, fill = "Condition"
  ) +
  theme_pub(10.5) +
  theme(legend.position = "top", plot.title.position = "plot")

paired_p <- suppressWarnings(wilcox.test(
  counts$potential_malignant_T_cell_percent[counts$condition == "Diagnosis"],
  counts$potential_malignant_T_cell_percent[counts$condition == "Relapse"],
  paired = TRUE,
  exact = TRUE,
  correct = FALSE
)$p.value)

p_prop <- ggplot(counts, aes(x = condition, y = potential_malignant_T_cell_percent, group = patient_pair)) +
  geom_line(color = "#9A9A9A", linewidth = 0.7) +
  geom_point(
    aes(color = condition, shape = condition, fill = condition),
    size = 4, stroke = 1.1
  ) +
  geom_text(data = counts %>% filter(condition == "Diagnosis"), aes(label = patient_pair), color = "#4A4A4A", size = 3.2, nudge_x = 0.08, nudge_y = 1.8, family = "Arial") +
  scale_color_manual(values = c("Diagnosis" = COL_DX, "Relapse" = COL_REL)) +
  scale_fill_manual(values = c("Diagnosis" = "white", "Relapse" = COL_REL)) +
  scale_shape_manual(values = c("Diagnosis" = 21, "Relapse" = 24)) +
  scale_y_continuous(limits = c(0, 78), breaks = seq(0, 75, 25), labels = label_percent(scale = 1)) +
  labs(
    title = "B  Paired cell proportions",
    subtitle = paste0("Paired Wilcoxon test, P = ", sprintf("%.2f", paired_p), "; N = 3 patients"),
    x = NULL, y = "Potential malignant T cells\n(% of candidate T-lineage cells)",
    color = "Condition", shape = "Condition", fill = "Condition"
  ) +
  theme_pub(10.5) +
  theme(legend.position = "none", plot.title.position = "plot")

fig1 <- p_count + p_prop + plot_layout(widths = c(1.65, 1))
save_plot(fig1, file.path(main_dir, "Fig01_cohort_cell_counts_and_paired_proportions.png"), 13.2, 5.2)

# -----------------------------------------------------------------------------
# Seurat-derived UMAP and paired STAT5 expression
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
})

obj <- readRDS(object_file)
potential_meta <- read_table("02_Full_Results", "potential_malignant_T_cells_metadata.csv")
primary_cells <- unique(potential_meta$cell)
umap <- as.data.frame(Embeddings(obj, reduction = "umap.harmony")) %>%
  rownames_to_column("cell")
colnames(umap)[2:3] <- c("UMAP_1", "UMAP_2")
meta <- obj@meta.data %>% rownames_to_column("cell")
umap_meta <- left_join(umap, meta, by = "cell") %>%
  mutate(
    population = case_when(
      cell %in% primary_cells ~ "Potential malignant T cells",
      broad_annotation == "T_lineage_candidate_malignant" ~ "Other candidate T-lineage cells",
      TRUE ~ "Other annotated cells"
    ),
    population = factor(population, levels = c("Other annotated cells", "Other candidate T-lineage cells", "Potential malignant T cells"))
  ) %>%
  arrange(population)

p_umap <- ggplot(umap_meta, aes(UMAP_1, UMAP_2, color = population)) +
  geom_point(size = 0.34, alpha = 0.78, stroke = 0) +
  scale_color_manual(values = c(
    "Other annotated cells" = "#D3D3D3",
    "Other candidate T-lineage cells" = "#E6AB02",
    "Potential malignant T cells" = "#8C1D18"
  )) +
  coord_equal() +
  labs(
    title = "Potential malignant T cells on the integrated UMAP",
    subtitle = "Primary population: CNV-high top 25% within candidate malignant T-lineage cells (3,207 cells)",
    x = "UMAP 1", y = "UMAP 2", color = NULL
  ) +
  theme_pub(11) +
  theme(legend.position = "right")
save_plot(p_umap, file.path(main_dir, "Fig02_primary_population_UMAP.png"), 8.4, 6.2)

# Retain the inferCNV panel after visual audit; rename it into the final figure series.
file.copy(
  file.path(old_figure_root, "01_Main_Figures", "infercnv_nonT_reference_final_heatmap.png"),
  file.path(main_dir, "Fig03_inferCNV_nonT_reference_heatmap.png"),
  overwrite = TRUE
)

# -----------------------------------------------------------------------------
# Figure 4. Focused KEGG GSEA
# -----------------------------------------------------------------------------
kegg_all <- read_table("02_Full_Results", "KEGG_GSEA_all_pathways_relapse_vs_diagnosis.csv")
kegg_targets <- c(
  "Ribosome biogenesis in eukaryotes",
  "Oxidative phosphorylation",
  "DNA replication",
  "Ribosome",
  "Spliceosome",
  "Alanine, aspartate and glutamate metabolism",
  "JAK-STAT signaling pathway",
  "HIF-1 signaling pathway",
  "Mitophagy - animal",
  "Proteasome",
  "Sphingolipid signaling pathway",
  "Autophagy - animal"
)
kegg_focus <- kegg_all %>%
  filter(pathway_name %in% kegg_targets, tested) %>%
  mutate(
    display = recode(pathway_name, "Mitophagy - animal" = "Mitophagy", "Autophagy - animal" = "Autophagy"),
    direction_label = ifelse(NES >= 0, "Relapse-up", "Relapse-down"),
    status = fdr_status(padj),
    q_label = paste0("q=", fmt_q(padj), " ", fdr_mark(padj)),
    display = factor(display, levels = rev(recode(kegg_targets, "Mitophagy - animal" = "Mitophagy", "Autophagy - animal" = "Autophagy")))
  )

p_kegg <- ggplot(kegg_focus, aes(y = display, x = NES)) +
  geom_segment(aes(x = 0, xend = NES, yend = display, color = direction_label), linewidth = 1.2, alpha = 0.55) +
  geom_point(aes(color = direction_label, shape = status), size = 4.2, stroke = 1.15) +
  geom_text(aes(label = q_label), hjust = ifelse(kegg_focus$NES >= 0, -0.12, 1.12), size = 3.25, family = "Arial") +
  geom_vline(xintercept = 0, color = "#777777", linewidth = 0.5) +
  scale_color_manual(values = c("Relapse-up" = COL_UP, "Relapse-down" = COL_DOWN)) +
  scale_shape_manual(values = status_shapes) +
  scale_x_continuous(limits = c(-1.55, 2.75), breaks = seq(-1.5, 2.5, 0.5)) +
  labs(
    title = "KEGG GSEA identifies coordinated relapse-associated programs",
    subtitle = "Paired patient-level pseudobulk model; * FDR < 0.05, # exploratory FDR < 0.10, ns not significant",
    x = "Normalized enrichment score (NES)", y = NULL,
    color = "Direction", shape = "Statistical evidence"
  ) +
  theme_pub(11) +
  theme(legend.position = "bottom", panel.grid.major.y = element_line(color = "#EFEFEF", linewidth = 0.35))
save_plot(p_kegg, file.path(main_dir, "Fig04_KEGG_GSEA_key_pathways.png"), 10.8, 7.2)

# -----------------------------------------------------------------------------
# Figure 5. Curated, nonredundant GO BP view
# -----------------------------------------------------------------------------
go_all <- read_table("02_Full_Results", "GO_BP_GSEA_all_terms_relapse_vs_diagnosis.csv") %>%
  filter(tested, !is.na(padj))
go_patterns <- tibble::tribble(
  ~display, ~regex,
  "Ribonucleoprotein complex biogenesis", "^ribonucleoprotein complex biogenesis$",
  "Cytoplasmic translation", "^cytoplasmic translation$",
  "Mitochondrial electron transport", "^mitochondrial ATP synthesis coupled electron transport$",
  "Aerobic electron transport chain", "^aerobic electron transport chain$",
  "DNA replication", "^DNA replication$|^DNA-dependent DNA replication$",
  "Cell-cycle DNA replication", "cell cycle DNA replication",
  "Spliceosomal mRNA processing", "^mRNA splicing, via spliceosome$|regulation of mRNA splicing, via spliceosome",
  "Proteasomal protein catabolism", "^proteasomal protein catabolic process$",
  "Type 2 mitophagy", "^type 2 mitophagy$",
  "Mitophagy", "^mitophagy$",
  "Autophagy", "^autophagy$",
  "Branched-chain amino acid catabolism", "^branched-chain amino acid catabolic process$"
)
go_focus <- bind_rows(lapply(seq_len(nrow(go_patterns)), function(i) {
  hit <- go_all %>% filter(grepl(go_patterns$regex[[i]], pathway_name, ignore.case = TRUE)) %>% arrange(padj) %>% slice_head(n = 1)
  if (nrow(hit) == 0) return(NULL)
  hit$display <- go_patterns$display[[i]]
  hit
})) %>%
  distinct(pathway_name, .keep_all = TRUE) %>%
  mutate(
    direction_label = ifelse(NES >= 0, "Relapse-up", "Relapse-down"),
    status = fdr_status(padj),
    q_label = paste0("q=", fmt_q(padj), " ", fdr_mark(padj)),
    display = factor(display, levels = rev(go_patterns$display[go_patterns$display %in% display]))
  )

p_go <- ggplot(go_focus, aes(y = display, x = NES)) +
  geom_segment(aes(x = 0, xend = NES, yend = display, color = direction_label), linewidth = 1.1, alpha = 0.52) +
  geom_point(aes(color = direction_label, shape = status), size = 4, stroke = 1.1) +
  geom_text(aes(label = q_label), hjust = ifelse(go_focus$NES >= 0, -0.12, 1.12), size = 3.15, family = "Arial") +
  geom_vline(xintercept = 0, color = "#777777", linewidth = 0.5) +
  scale_color_manual(values = c("Relapse-up" = COL_UP, "Relapse-down" = COL_DOWN)) +
  scale_shape_manual(values = status_shapes) +
  scale_x_continuous(limits = c(-1.9, 2.65), breaks = seq(-1.5, 2.5, 0.5)) +
  labs(
    title = "GO biological processes support translation and respiratory remodeling",
    subtitle = "Selected nonredundant terms; the complete 4,746 tested terms remain available as CSV",
    x = "Normalized enrichment score (NES)", y = NULL,
    color = "Direction", shape = "Statistical evidence"
  ) +
  theme_pub(11) +
  theme(legend.position = "bottom", panel.grid.major.y = element_line(color = "#EFEFEF", linewidth = 0.35))
save_plot(p_go, file.path(main_dir, "Fig05_GO_BP_GSEA_curated_programs.png"), 10.9, 7.1)

# -----------------------------------------------------------------------------
# Figure 6. Leading-edge heatmap with the correct pathway-level statistics
# -----------------------------------------------------------------------------
heat_genes <- read_table("02_Full_Results", "proteasome_mitophagy_GSEA_heatmap_selected_genes.csv") %>%
  group_by(pathway_name) %>% arrange(desc(rank_stat), .by_group = TRUE) %>% slice_head(n = 8) %>% ungroup() %>%
  mutate(
    pathway_short = recode(pathway_name, "Mitophagy - animal" = "Mitophagy"),
    gene_order = row_number()
  )
heat_expr <- read_table("02_Full_Results", "proteasome_mitophagy_sample_pseudobulk_logCPM.csv") %>%
  filter(gene %in% heat_genes$gene) %>%
  mutate(
    condition = factor(condition_simple, levels = c("diagnosis", "relapse"), labels = c("Diagnosis", "Relapse")),
    sample_pub = paste0(patient_pair, "\n", ifelse(condition == "Diagnosis", "Dx", "Rel")),
    sample_pub = factor(sample_pub, levels = c("P1\nDx", "P1\nRel", "P2\nDx", "P2\nRel", "P3\nDx", "P3\nRel"))
  ) %>%
  group_by(gene) %>%
  mutate(z = as.numeric(scale(pseudobulk_logCPM))) %>%
  ungroup() %>%
  left_join(heat_genes %>% select(gene, pathway_short, gene_order), by = "gene") %>%
  mutate(gene = factor(gene, levels = rev(unique(heat_genes$gene))))

cond_ann <- distinct(heat_expr, sample_pub, condition)
p_cond <- ggplot(cond_ann, aes(sample_pub, 1, fill = condition)) +
  geom_tile(color = "white", linewidth = 0.8) +
  scale_fill_manual(values = c("Diagnosis" = COL_DX, "Relapse" = COL_REL)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, fill = "Condition") +
  theme_void(base_family = "Arial") +
  theme(legend.position = "top", plot.margin = margin(0, 5, 0, 5))

p_heat <- ggplot(heat_expr, aes(sample_pub, gene, fill = z)) +
  geom_tile(color = "white", linewidth = 0.45) +
  facet_grid(pathway_short ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-2, 2), oob = squish) +
  labs(x = NULL, y = NULL, fill = "Row z-score") +
  theme_minimal(base_size = 10.5, base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold", color = "#202020"),
    axis.text.y = element_text(face = "italic", color = "#202020"),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold", margin = margin(r = 8)),
    strip.background = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(0, 5, 4, 5)
  )

path_summary <- heat_genes %>%
  distinct(pathway_short, pathway_NES, pathway_FDR) %>%
  mutate(
    status = fdr_status(pathway_FDR),
    mark = fdr_mark(pathway_FDR),
    label = paste0("NES ", sprintf("%.2f", pathway_NES), "\nq=", fmt_q(pathway_FDR), "  ", mark),
    pathway_short = factor(pathway_short, levels = c("Proteasome", "Mitophagy"))
  )
p_path <- ggplot(path_summary, aes(pathway_NES, pathway_short)) +
  geom_segment(aes(x = 0, xend = pathway_NES, yend = pathway_short), color = "#B9B9B9", linewidth = 1.1) +
  geom_point(aes(shape = status), color = COL_REL, size = 4.2, stroke = 1.1) +
  geom_text(aes(label = label), nudge_x = 0.08, hjust = 0, size = 3.35, family = "Arial") +
  scale_shape_manual(values = status_shapes) +
  scale_x_continuous(limits = c(0, 2.35), breaks = c(0, 1, 2)) +
  labs(title = "Pathway GSEA", x = "NES", y = NULL, shape = "Evidence") +
  theme_pub(10.5) +
  theme(legend.position = "bottom", axis.text.y = element_text(face = "bold"))

fig6_left <- p_cond / p_heat + plot_layout(heights = c(0.55, 8.5), guides = "collect")
fig6 <- fig6_left | p_path
fig6 <- fig6 + plot_annotation(
  title = "Mitophagy and proteasome leading-edge expression across paired samples",
  subtitle = "No selected gene passed gene-level FDR < 0.05; symbols report pathway-level evidence (* FDR < 0.05, # exploratory FDR < 0.10, ns not significant)",
  theme = theme(plot.title = element_text(face = "bold", size = 15, family = "Arial"), plot.subtitle = element_text(size = 10.5, family = "Arial", color = "#3F3F3F"))
)
save_plot(fig6, file.path(main_dir, "Fig06_mitophagy_proteasome_leading_edge_heatmap.png"), 12.8, 7.8)

# -----------------------------------------------------------------------------
# Figure 7. Cross-dataset pathway comparison
# -----------------------------------------------------------------------------
cross <- read_table("01_Key_Results", "cross_dataset_STAT5B_N642H_vs_relapse_GSEA_key_pathways.csv") %>%
  mutate(
    dataset_short = recode(dataset,
      "Relapse scRNA-seq" = "Paired relapse scRNA-seq",
      "STAT5B N642H bulk RNA-seq" = "STAT5B N642H bulk RNA-seq"
    ),
    status = fdr_status(FDR),
    q_label = paste0("q=", fmt_q(FDR), " ", fdr_mark(FDR)),
    pathway = factor(pathway, levels = rev(unique(pathway)))
  )
p_cross <- ggplot(cross, aes(NES, pathway)) +
  geom_vline(xintercept = 0, color = "#777777", linewidth = 0.45) +
  geom_segment(aes(x = 0, xend = NES, yend = pathway, color = dataset_short), linewidth = 1.05, alpha = 0.5) +
  geom_point(aes(color = dataset_short, shape = status), size = 4, stroke = 1.15) +
  geom_text(aes(label = q_label), hjust = ifelse(cross$NES >= 0, -0.12, 1.12), size = 2.95, family = "Arial") +
  facet_wrap(~dataset_short, nrow = 1) +
  scale_color_manual(values = c("Paired relapse scRNA-seq" = COL_DX, "STAT5B N642H bulk RNA-seq" = COL_BULK)) +
  scale_shape_manual(values = status_shapes) +
  scale_x_continuous(limits = c(-2.55, 3.15), breaks = seq(-2, 3, 1)) +
  labs(
    title = "Pathway concordance with constitutively active STAT5B",
    subtitle = "Same pathway direction supports association, not direct causality; * FDR < 0.05, # exploratory FDR < 0.10",
    x = "Normalized enrichment score (NES)", y = NULL, shape = "Statistical evidence", color = "Dataset"
  ) +
  theme_pub(10.5) +
  theme(legend.position = "bottom", panel.grid.major.y = element_line(color = "#EFEFEF", linewidth = 0.35))
save_plot(p_cross, file.path(main_dir, "Fig07_cross_dataset_GSEA_concordance.png"), 13, 6.8)

# -----------------------------------------------------------------------------
# Figure 8. Same-direction KEGG overlap at FDR < 0.10
# -----------------------------------------------------------------------------
venn_counts <- read_table("01_Key_Results", "cross_dataset_KEGG_FDR0.10_same_direction_venn_counts.csv")
circle_points <- function(cx, cy, r, set_name) {
  th <- seq(0, 2 * pi, length.out = 300)
  data.frame(x = cx + r * cos(th), y = cy + r * sin(th), set_name = set_name)
}
make_venn <- function(direction_name, row_data) {
  circles <- bind_rows(circle_points(0.85, 1, 0.78, "scRNA-seq"), circle_points(1.75, 1, 0.78, "STAT5B bulk"))
  ggplot(circles, aes(x, y, group = set_name, fill = set_name, color = set_name)) +
    geom_polygon(alpha = 0.34, linewidth = 1.1) +
    annotate("text", x = 0.35, y = 1, label = row_data$sc_only, size = 6, family = "Arial", fontface = "bold") +
    annotate("text", x = 1.30, y = 1, label = row_data$shared, size = 6.2, family = "Arial", fontface = "bold") +
    annotate("text", x = 2.25, y = 1, label = row_data$bulk_only, size = 6, family = "Arial", fontface = "bold") +
    annotate("text", x = 0.55, y = 1.83, label = "Relapse\nscRNA-seq", size = 3.5, family = "Arial", color = COL_DX, fontface = "bold") +
    annotate("text", x = 2.05, y = 1.83, label = "STAT5B N642H\nbulk RNA-seq", size = 3.5, family = "Arial", color = "#9C6B00", fontface = "bold") +
    coord_equal(xlim = c(-0.05, 2.65), ylim = c(0.05, 2.15), clip = "off") +
    scale_fill_manual(values = c("scRNA-seq" = COL_DX, "STAT5B bulk" = COL_BULK)) +
    scale_color_manual(values = c("scRNA-seq" = COL_DX, "STAT5B bulk" = COL_BULK)) +
    labs(title = paste0(direction_name, " pathways")) +
    theme_void(base_family = "Arial") +
    theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5), legend.position = "none", plot.margin = margin(8, 8, 8, 8))
}
p_up <- make_venn("Up-enriched", venn_counts %>% filter(direction == "up"))
p_down <- make_venn("Down-enriched", venn_counts %>% filter(direction == "down"))
fig8 <- p_up + p_down + plot_annotation(
  title = "Same-direction KEGG overlap at exploratory FDR < 0.10",
  subtitle = "Twenty-eight pathways were shared in the upward direction; none were shared in the downward direction",
  theme = theme(plot.title = element_text(face = "bold", size = 15, family = "Arial"), plot.subtitle = element_text(size = 10.5, family = "Arial", color = "#3F3F3F"))
)
save_plot(fig8, file.path(main_dir, "Fig08_cross_dataset_KEGG_overlap_venn.png"), 11.6, 5.7)

# -----------------------------------------------------------------------------
# Supplementary annotation UMAP and audited marker dot plot
# -----------------------------------------------------------------------------
ann_palette <- c(
  "B-cell-like_contaminant" = "#0072B2",
  "Quiescent_CD44_high_TALL_like" = "#CC79A7",
  "Cycling_TALL_blast_like" = "#D55E00",
  "Cytotoxic_T_or_NK_like" = "#009E73",
  "ER_stress_TALL_like" = "#AA4499",
  "Erythroid-like_contaminant" = "#E64B35",
  "Immature_HSC_like_TALL" = "#56B4E9",
  "Myeloid-like_contaminant" = "#999933",
  "pDC-like_contaminant" = "#E69F00"
)
ann_levels <- names(ann_palette)
ann_labels <- c(
  "B cell-like", "CD44-high quiescent T-ALL-like", "Cycling T-ALL-like", "Cytotoxic T/NK-like",
  "ER-stress T-ALL-like", "Erythroid-like", "Immature/HSC-like T-ALL", "Myeloid-like", "pDC-like"
)
p_ann <- umap_meta %>%
  mutate(major_annotation = factor(major_annotation, levels = ann_levels)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = major_annotation)) +
  geom_point(size = 0.32, alpha = 0.78, stroke = 0) +
  scale_color_manual(values = ann_palette, labels = ann_labels, drop = FALSE) +
  coord_equal() +
  labs(
    title = "Marker-based major cell annotation",
    subtitle = "Distinct colors are used only for biologically separate annotation classes",
    x = "UMAP 1", y = "UMAP 2", color = "Annotation"
  ) +
  theme_pub(10.5) +
  theme(legend.position = "right")
save_plot(p_ann, file.path(supp_dir, "FigS01_marker_based_annotation_UMAP.png"), 9.2, 6.5)

file.copy(
  file.path(old_figure_root, "01_Main_Figures", "annotation_lineage_marker_dotplot_publication.png"),
  file.path(supp_dir, "FigS02_lineage_marker_dotplot.png"),
  overwrite = TRUE
)

# -----------------------------------------------------------------------------
# Supplementary STAT5 figure: paired expression plus downstream GSEA
# -----------------------------------------------------------------------------
counts_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")
primary_cells <- intersect(primary_cells, colnames(counts_mat))
sample_meta <- potential_meta %>% distinct(sample_id, patient_pair, condition_simple)
sample_ids <- sample_meta$sample_id
pseudobulk <- sapply(sample_ids, function(sid) {
  cs <- intersect(potential_meta$cell[potential_meta$sample_id == sid], colnames(counts_mat))
  Matrix::rowSums(counts_mat[, cs, drop = FALSE])
})
colnames(pseudobulk) <- sample_ids
dge <- DGEList(pseudobulk)
dge <- calcNormFactors(dge)
log_cpm <- cpm(dge, log = TRUE, prior.count = 2)
stat5_expr <- as.data.frame(log_cpm[c("STAT5A", "STAT5B"), , drop = FALSE]) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "logCPM") %>%
  left_join(sample_meta, by = "sample_id") %>%
  mutate(condition = factor(condition_simple, levels = c("diagnosis", "relapse"), labels = c("Diagnosis", "Relapse")))
stat5_deg <- read_table("01_Key_Results", "STAT5A_STAT5B_pseudobulk_DEG.csv") %>%
  mutate(label = paste0("log2FC ", sprintf("%+.2f", logFC), "\nnominal P=", fmt_q(PValue), ", FDR=", fmt_q(FDR)))

p_stat5_expr <- ggplot(stat5_expr, aes(condition, logCPM, group = patient_pair)) +
  geom_line(color = "#9A9A9A", linewidth = 0.65) +
  geom_point(aes(color = condition, shape = condition, fill = condition), size = 3.5, stroke = 1) +
  facet_wrap(~gene, scales = "free_y") +
  scale_color_manual(values = c("Diagnosis" = COL_DX, "Relapse" = COL_REL)) +
  scale_fill_manual(values = c("Diagnosis" = "white", "Relapse" = COL_REL)) +
  scale_shape_manual(values = c("Diagnosis" = 21, "Relapse" = 24)) +
  labs(title = "A  Paired STAT5 expression", subtitle = "Patient-level pseudobulk logCPM", x = NULL, y = "log2 counts per million", color = NULL, shape = NULL, fill = NULL) +
  theme_pub(10.2) +
  theme(legend.position = "top")

stat5_sets <- read_table("01_Key_Results", "STAT5_downstream_gene_set_GSEA_relapse_vs_diagnosis.csv") %>%
  mutate(
    display = recode(pathway,
      "STAT5_canonical_feedback_survival" = "Canonical feedback / survival",
      "STAT5_T_cell_cytokine_program" = "T-cell cytokine program",
      "STAT5_leukemia_relevant_core" = "Leukemia-relevant core"
    ),
    status = fdr_status(padj),
    q_label = paste0("q=", fmt_q(padj), " ", fdr_mark(padj)),
    display = factor(display, levels = rev(display))
  )
p_stat5_sets <- ggplot(stat5_sets, aes(NES, display)) +
  geom_segment(aes(x = 0, xend = NES, yend = display), color = "#BDBDBD", linewidth = 1) +
  geom_point(aes(shape = status), size = 4, color = COL_REL, stroke = 1.1) +
  geom_text(aes(label = q_label), nudge_x = 0.06, hjust = 0, size = 3.2, family = "Arial") +
  scale_shape_manual(values = status_shapes) +
  scale_x_continuous(limits = c(0, 1.72), breaks = c(0, 0.5, 1, 1.5)) +
  labs(title = "B  STAT5 downstream gene-set GSEA", subtitle = "All tested sets were not significant", x = "NES", y = NULL, shape = "Evidence") +
  theme_pub(10.2) +
  theme(legend.position = "bottom")

stat5_note <- ggplot() +
  annotate("text", x = 0, y = 1, label = paste0(stat5_deg$gene, ":  ", stat5_deg$label, collapse = "\n\n"), hjust = 0, vjust = 1, family = "Arial", size = 3.7, lineheight = 1.15) +
  xlim(0, 1) + ylim(0, 1) + theme_void() + labs(title = "Gene-level paired model") +
  theme(plot.title = element_text(face = "bold", family = "Arial", size = 11.5))

fig_stat5 <- p_stat5_expr | (p_stat5_sets / stat5_note + plot_layout(heights = c(1.3, 1)))
fig_stat5 <- fig_stat5 + plot_annotation(
  title = "STAT5A expression increased, but STAT5 pathway evidence is not FDR-significant",
  subtitle = "Nominal P values and FDR are shown together to avoid overstating a two-gene view",
  theme = theme(plot.title = element_text(face = "bold", size = 14, family = "Arial"), plot.subtitle = element_text(size = 10, family = "Arial", color = "#3F3F3F"))
)
save_plot(fig_stat5, file.path(supp_dir, "FigS03_STAT5_paired_expression_and_downstream_GSEA.png"), 12.2, 6.8)

# -----------------------------------------------------------------------------
# Supplementary custom mitochondrial gene sets
# -----------------------------------------------------------------------------
custom_mito <- read_table("04_Method_Inputs", "custom_mitochondrial_gene_set_GSEA.csv") %>%
  mutate(
    display = recode(pathway,
      "mitochondrial_oxphos_core" = "Mitochondrial OXPHOS core",
      "mitochondrial_stress_dynamics" = "Mitochondrial stress / dynamics",
      "mitophagy_core" = "Mitophagy core",
      "proteasome_core" = "Proteasome core",
      "mitochondrial_translation" = "Mitochondrial translation",
      "autophagy_core" = "Autophagy core"
    ),
    direction_label = ifelse(NES >= 0, "Relapse-up", "Relapse-down"),
    status = fdr_status(padj),
    q_label = paste0("q=", fmt_q(padj), " ", fdr_mark(padj)),
    display = factor(display, levels = rev(display))
  )
p_custom <- ggplot(custom_mito, aes(NES, display)) +
  geom_vline(xintercept = 0, color = "#777777", linewidth = 0.45) +
  geom_segment(aes(x = 0, xend = NES, yend = display, color = direction_label), linewidth = 1.1, alpha = 0.55) +
  geom_point(aes(color = direction_label, shape = status), size = 4.2, stroke = 1.1) +
  geom_text(aes(label = q_label), hjust = ifelse(custom_mito$NES >= 0, -0.12, 1.12), size = 3.2, family = "Arial") +
  scale_color_manual(values = c("Relapse-up" = COL_UP, "Relapse-down" = COL_DOWN)) +
  scale_shape_manual(values = status_shapes) +
  scale_x_continuous(limits = c(-1.55, 2.8), breaks = seq(-1.5, 2.5, 0.5)) +
  labs(
    title = "Curated mitochondrial and proteostasis gene-set GSEA",
    subtitle = "Custom core sets are supportive analyses and are distinct from the KEGG gene sets",
    x = "Normalized enrichment score (NES)", y = NULL, color = "Direction", shape = "Statistical evidence"
  ) +
  theme_pub(11) + theme(legend.position = "bottom")
save_plot(p_custom, file.path(supp_dir, "FigS04_custom_mitochondrial_gene_set_GSEA.png"), 9.8, 5.6)

# -----------------------------------------------------------------------------
# Supplementary sphingolipid paired-direction matrix
# -----------------------------------------------------------------------------
sph <- read_table("01_Key_Results", "sphingolipid_leading_edge_gene_paired_directions.csv") %>%
  arrange(desc(abs(logFC))) %>% slice_head(n = 16) %>%
  select(gene, P1_delta, P2_delta, P3_delta, paired_direction_consistency, FDR) %>%
  pivot_longer(matches("^P[1-3]_delta$"), names_to = "patient", values_to = "paired_delta") %>%
  mutate(patient = sub("_delta", "", patient), gene = factor(gene, levels = rev(unique(gene))))
sph_summary <- read_table("01_Key_Results", "sphingolipid_KEGG_GSEA_summary.csv")
sph_text <- sph_summary %>% filter(pathway_name %in% c("Sphingolipid signaling pathway", "Sphingolipid metabolism")) %>%
  transmute(txt = paste0(pathway_name, ": NES ", sprintf("%.2f", NES), ", q=", fmt_q(GSEA_FDR), " ns")) %>% pull(txt) %>% paste(collapse = "\n")
p_sph <- ggplot(sph, aes(patient, gene, fill = paired_delta)) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(aes(label = sprintf("%+.1f", paired_delta)), size = 3, family = "Arial") +
  scale_fill_gradient2(low = COL_DOWN, mid = "white", high = COL_UP, midpoint = 0, limits = c(-3, 3), oob = squish) +
  labs(
    title = "Sphingolipid leading-edge genes show mixed but often concordant paired directions",
    subtitle = paste0(sph_text, "\nNo displayed gene passed gene-level FDR < 0.05."),
    x = "Patient pair", y = NULL, fill = "Relapse minus\ndiagnosis logCPM"
  ) +
  theme_pub(10.5) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(face = "italic"), legend.position = "right")
save_plot(p_sph, file.path(supp_dir, "FigS05_sphingolipid_paired_direction_heatmap.png"), 9.3, 7.2)

# -----------------------------------------------------------------------------
# Sensitivity heatmaps. Primary top-25% rows are outlined.
# -----------------------------------------------------------------------------
compartment_labels <- c(
  "CNV_high_top25_candidate" = "Primary: CNV-high top 25%",
  "CNV_high_top50_candidate" = "CNV-high top 50%",
  "broad_T_lineage_candidate_malignant" = "All candidate T-lineage",
  "Cycling_TALL_blast_like" = "Cycling T-ALL-like",
  "Quiescent_CD44_high_TALL_like" = "CD44-high quiescent T-ALL-like"
)
sens_heatmap <- function(dat, pathway_col, pathway_levels, filename, title, width = 10.5, height = 6.2) {
  dat2 <- dat %>%
    filter(compartment %in% names(compartment_labels), .data[[pathway_col]] %in% names(pathway_levels), tested, !is.na(NES)) %>%
    mutate(
      compartment_display = factor(recode(compartment, !!!compartment_labels), levels = unname(compartment_labels)),
      pathway_display = factor(recode(.data[[pathway_col]], !!!pathway_levels), levels = unname(pathway_levels)),
      cell_label = paste0(sprintf("%.2f", NES), " ", fdr_mark(padj)),
      is_primary = compartment == "CNV_high_top25_candidate"
    )
  p <- ggplot(dat2, aes(pathway_display, compartment_display, fill = NES)) +
    geom_tile(color = "white", linewidth = 0.7) +
    geom_tile(data = dat2 %>% filter(is_primary), fill = NA, color = "#1A1A1A", linewidth = 1.05) +
    geom_text(aes(label = cell_label), size = 3.15, family = "Arial") +
    scale_fill_gradient2(low = COL_DOWN, mid = "white", high = COL_UP, midpoint = 0, limits = c(-2.4, 2.4), oob = squish) +
    labs(
      title = title,
      subtitle = "Cells show NES; * FDR < 0.05, # exploratory FDR < 0.10, ns not significant. Black outline marks the primary population.",
      x = NULL, y = NULL, fill = "NES"
    ) +
    theme_pub(10.5) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      axis.text.y = element_text(face = ifelse(levels(dat2$compartment_display) == "Primary: CNV-high top 25%", "bold", "plain")),
      legend.position = "right"
    )
  save_plot(p, filename, width, height)
}

sens_custom <- read_table("03_Sensitivity_Checks", "sensitivity_readout_custom_gsea_valid_compartments.csv")
sens_heatmap(
  sens_custom, "pathway_name",
  c(
    "mitochondrial_oxphos_core" = "OXPHOS core",
    "mitochondrial_translation" = "Mitochondrial translation",
    "mitochondrial_stress_dynamics" = "Stress / dynamics",
    "mitophagy_core" = "Mitophagy core",
    "proteasome_core" = "Proteasome core",
    "autophagy_core" = "Autophagy core"
  ),
  file.path(sens_dir, "FigSens01_custom_gene_sets_across_cell_definitions.png"),
  "Sensitivity analysis of curated mitochondrial and proteostasis gene sets"
)

sens_kegg <- read_table("03_Sensitivity_Checks", "sensitivity_readout_key_kegg_focus.csv")
sens_heatmap(
  sens_kegg, "pathway_name",
  c(
    "Oxidative phosphorylation" = "OXPHOS",
    "Ribosome" = "Ribosome",
    "DNA replication" = "DNA replication",
    "JAK-STAT signaling pathway" = "JAK-STAT",
    "Mitophagy - animal" = "Mitophagy",
    "Proteasome" = "Proteasome",
    "Autophagy - animal" = "Autophagy"
  ),
  file.path(sens_dir, "FigSens02_KEGG_across_cell_definitions.png"),
  "Sensitivity analysis of key KEGG pathways",
  11.2, 6.3
)

sens_go <- read_table("03_Sensitivity_Checks", "sensitivity_readout_go_mitophagy_rows.csv")
sens_heatmap(
  sens_go, "pathway_name",
  c(
    "type 2 mitophagy" = "Type 2 mitophagy",
    "positive regulation of mitophagy" = "Positive regulation of mitophagy",
    "regulation of mitophagy" = "Regulation of mitophagy",
    "mitophagy" = "Mitophagy"
  ),
  file.path(sens_dir, "FigSens03_GO_mitophagy_across_cell_definitions.png"),
  "GO mitophagy terms are not independently significant across cell definitions",
  10.5, 5.8
)

# -----------------------------------------------------------------------------
# Final indexes and figure-audit record
# -----------------------------------------------------------------------------
figure_index <- tribble(
  ~priority, ~file, ~what_it_shows, ~statistical_message,
  "Main", "01_Main_Figures/Fig01_cohort_cell_counts_and_paired_proportions.png", "Counts and paired proportions of the primary population", "Paired Wilcoxon P shown; N=3 patients",
  "Main", "01_Main_Figures/Fig02_primary_population_UMAP.png", "Only the primary CNV-high top-25% definition is highlighted", "Cell-definition visualization; not an outcome test",
  "Main", "01_Main_Figures/Fig03_inferCNV_nonT_reference_heatmap.png", "inferCNV pattern relative to non-T references", "Basis for the CNV burden score",
  "Main", "01_Main_Figures/Fig04_KEGG_GSEA_key_pathways.png", "Focused paired-pseudobulk KEGG GSEA", "FDR evidence encoded by both shape and text",
  "Main", "01_Main_Figures/Fig05_GO_BP_GSEA_curated_programs.png", "Curated nonredundant GO BP programs", "Full tested table retained; plot avoids redundant GO terms",
  "Main", "01_Main_Figures/Fig06_mitophagy_proteasome_leading_edge_heatmap.png", "Paired sample expression of selected leading-edge genes", "Gene-level FDR absent; pathway-level FDR shown separately",
  "Main", "01_Main_Figures/Fig07_cross_dataset_GSEA_concordance.png", "Relapse scRNA-seq compared with STAT5B N642H bulk RNA-seq", "Association only; no causal claim",
  "Main", "01_Main_Figures/Fig08_cross_dataset_KEGG_overlap_venn.png", "Same-direction KEGG overlap at FDR < 0.10", "28 shared up pathways and 0 shared down pathways",
  "Supplementary", "02_Supplementary_Figures/FigS01_marker_based_annotation_UMAP.png", "Major cell annotation", "Colorblind-aware categorical palette",
  "Supplementary", "02_Supplementary_Figures/FigS02_lineage_marker_dotplot.png", "Marker support for annotations", "Expression and detection are encoded separately",
  "Supplementary", "02_Supplementary_Figures/FigS03_STAT5_paired_expression_and_downstream_GSEA.png", "STAT5A/B paired expression and downstream GSEA", "Nominal P and FDR are displayed together",
  "Supplementary", "02_Supplementary_Figures/FigS04_custom_mitochondrial_gene_set_GSEA.png", "Curated mitochondrial and proteostasis sets", "Custom-set evidence separated from KEGG evidence",
  "Supplementary", "02_Supplementary_Figures/FigS05_sphingolipid_paired_direction_heatmap.png", "Per-patient sphingolipid leading-edge directions", "Pathways and genes are not FDR-significant",
  "Sensitivity", "03_Sensitivity_Checks/FigSens01_custom_gene_sets_across_cell_definitions.png", "Custom gene-set robustness", "Primary population outlined; alternatives labeled sensitivity only",
  "Sensitivity", "03_Sensitivity_Checks/FigSens02_KEGG_across_cell_definitions.png", "KEGG robustness", "Primary population outlined; alternatives labeled sensitivity only",
  "Sensitivity", "03_Sensitivity_Checks/FigSens03_GO_mitophagy_across_cell_definitions.png", "GO mitophagy robustness", "GO mitophagy terms remain non-significant"
)
write.csv(figure_index, file.path(output_root, "FIGURE_INDEX.csv"), row.names = FALSE, na = "")

audit <- tribble(
  ~issue_in_previous_package, ~correction_in_publication_suite,
  "Diagnosis and relapse were sometimes distinguished by color alone", "Condition is distinguished by color plus point shape and fill",
  "Nominal P values could be mistaken for FDR significance", "Nominal P and FDR are printed together; pathway and gene statistics are separated",
  "GO figures contained many redundant or biologically distracting terms", "Main GO figure uses a documented nonredundant selection; complete results remain in CSV",
  "The primary UMAP also displayed obsolete top-50% thresholds", "Main UMAP highlights only the pre-specified CNV-high top-25% population",
  "Proteasome and mitophagy heatmaps mixed pathway and gene-level significance", "No gene-level stars are shown because no selected gene passes FDR; pathway FDR has its own panel",
  "Sensitivity plots used underscored analysis labels and similar encodings", "Human-readable labels, a common NES scale and a black outline for the primary population are used",
  "Several main figures duplicated the same GSEA result", "Each biological conclusion now has one main figure; redundant and thresholded ORA plots are removed"
)
write.csv(audit, file.path(output_root, "FIGURE_AUDIT.csv"), row.names = FALSE, na = "")

message("Publication figure suite written to: ", output_root)
