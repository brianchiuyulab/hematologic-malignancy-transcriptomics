# Code availability

This folder contains the analysis code used for GSE262271 and the external STAT5B N642H comparison. Scripts are numbered by conceptual run order. Scripts 01 to 15 preserve the full analysis provenance. Script 16 is the canonical publication-figure renderer and is the only script that should write the final figure suite.

## Statistical design

The primary contrast is relapse versus diagnosis within the same patient. Cells were aggregated to patient-sample pseudobulk counts. The edgeR design included both patient and condition:

~ patient_pair + condition_simple

The condition coefficient estimates the average within-patient relapse effect. GSEA used a genome-wide ranking defined as:

sign(logFC) × sqrt(F)

This is not an independent-cell t test and cells are not treated as biological replicates.

## Primary population

All main disease-comparison scripts use the CNV-high top 25% subset within candidate malignant T-lineage cells. Broader definitions occur only in script 12 and are labeled sensitivity analyses.

## Script map

| Step | Script | Purpose |
|---|---|---|
| 1 | 01_full_preprocessing_and_seurat.R | Raw import, QC, integration and Seurat object |
| 2 | 02_downstream_annotation_and_infercnv_inputs.R | Annotation, pseudobulk foundations and inferCNV inputs |
| 3 | 03_run_infercnv_nonT_reference.R | inferCNV with non-T reference groups |
| 4 | 04_summarize_infercnv_burden.R | Cell-level CNV burden and CNV-high subsets |
| 5 | 05_annotation_QC_figures.R | Annotation and CNV definition figures |
| 6 | 06_potential_malignant_T_KEGG_ORA_GSEA.R | Primary KEGG ORA and GSEA |
| 7 | 07_potential_malignant_T_GO_BP_GSEA.R | GO biological-process GSEA |
| 8 | 08_potential_malignant_T_STAT5_GSEA.R | Custom STAT5 downstream gene sets |
| 9 | 09_mitophagy_autophagy_paired_consistency.R | Per-patient leading-edge direction audit |
| 10 | 10_proteasome_mitophagy_publication_figure.R | Publication-style heatmap and paired scores |
| 11 | 11_sphingolipid_paired_analysis.R | Sphingolipid-focused analysis |
| 12 | 12_sensitivity_analysis.R | Broader-cell-definition robustness checks |
| 13–15 | cross-dataset scripts | GSE218858 STAT5B N642H validation and overlap |
| 16 | 16_publication_figure_suite.R | Canonical, colorblind-aware final figure suite and figure audit |

The editable handoff deck is built by `presentation/build_handoff_ppt.mjs`. It uses the final PNG figures, creates native PowerPoint tables for the study design and result summaries, and is validated for slide count, geometry and font consistency before release. Rebuilding it requires the OpenAI Artifact Tool presentation runtime; the reviewed `.pptx` is included for users who do not have that runtime.

## Reproducing figures without altering the clean package

1. Run `run_final_figures.ps1` to regenerate the publication figure suite from the frozen result tables and final Seurat object.
2. The wrapper writes to an ignored `_reproduction_run/publication_figures` folder, so it cannot overwrite the reviewed final figures.
3. For a complete raw-to-result rerun, first run `prepare_reproduction_workspace.ps1`, then use scripts 01 to 12 inside `_reproduction_run`. inferCNV requires JAGS and the R package infercnv.
4. Cross-dataset scripts 13 to 15 require the separately analyzed GSE218858 folder. Their frozen outputs are already included in the compact tables used by script 16.

Historical scripts contain Windows defaults from the original workstation. The wrappers always pass explicit paths, which should be preferred over those defaults.

## Required R packages

Core packages include Seurat, SeuratObject, Matrix, harmony, edgeR, fgsea, dplyr, tidyr, tibble, ggplot2, patchwork, AnnotationDbi, org.Hs.eg.db, GO.db and infercnv. JAGS is required for the inferCNV installation used here. Exact session records are under environment.

## GitHub availability

The hematologic malignancy project repository contains the code, README, final figures, compact CSV tables and handoff PPT. `03_Raw_Data` and `04_R_Objects` are excluded. Reproduction therefore starts from the frozen compact tables unless the user supplies the GEO raw archive and final R objects locally.
