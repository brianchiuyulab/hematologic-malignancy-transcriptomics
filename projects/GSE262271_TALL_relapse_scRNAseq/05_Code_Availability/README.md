# Code availability

This folder contains the analysis code used for GSE262271 and the external STAT5B N642H comparison. Scripts are numbered by conceptual run order. Scripts 01 to 15 preserve the full analysis provenance. Script 16 is the canonical publication-figure renderer. Scripts 17 to 20 write separate sensitivity, fixed-secondary or publication-compatibility outputs without changing the frozen primary result.

## Statistical design

The primary contrast is relapse versus diagnosis within the same patient. Cells were aggregated to patient-sample pseudobulk counts. The edgeR design included both patient and condition:

~ patient_pair + condition_simple

The condition coefficient estimates the average within-patient relapse effect. GSEA used a genome-wide ranking defined as:

sign(logFC) × sqrt(F)

This is not an independent-cell t test and cells are not treated as biological replicates.

## Primary population

All main disease-comparison scripts use the CNV-high top 25% subset within candidate malignant T-lineage cells. Broader definitions or alternative parameters occur only in scripts 12 and 17 and are labeled sensitivity analyses.

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
| 17 | 17_robustness_sensitivity_audit.R | Predefined CNV, QC, filter, model, ranking and leave-one-patient-out robustness grid |
| 18 | 18_targeted_GSEA_multiverse_sensitivity.R | Exhaustive targeted grid for official KEGG Proteasome and Mitophagy across 528 population/QC/model/ranking specifications |
| 19 | 19_fixed_secondary_within_sample_top25_KEGG.R | Reader-facing complete-KEGG tables and mechanistic dot plot for source-paper QC plus CNV-high top 25% within each sample |
| 20 | 20_publication_compatible_candidate_full_KEGG.R | Full-KEGG audit of 10 publication-defensible population/QC candidates, with global FDR reported for Proteasome and Mitophagy |

The editable handoff deck is built by `presentation/build_handoff_ppt.mjs`. It uses the final PNG figures, creates native PowerPoint tables for the study design and result summaries, and is validated for slide count, geometry and font consistency before release. Rebuilding it requires the OpenAI Artifact Tool presentation runtime; the reviewed `.pptx` is included for users who do not have that runtime.

## Reproducing figures without altering the clean package

1. Run `run_final_figures.ps1` to regenerate the publication figure suite from the frozen result tables and final Seurat object.
2. The wrapper writes to an ignored `_reproduction_run/publication_figures` folder, so it cannot overwrite the reviewed final figures.
3. For a complete raw-to-result rerun, first run `prepare_reproduction_workspace.ps1`, then use scripts 01 to 12 inside `_reproduction_run`. inferCNV requires JAGS and the R package infercnv.
4. Cross-dataset scripts 13 to 15 require the separately analyzed GSE218858 folder. Their frozen outputs are already included in the compact tables used by script 16.
5. Run script 17 from the project root to regenerate the sensitivity audit. Use `--reuse-results true` only to redraw its figures from an already completed grid.
6. Run script 18 with `--nperm 5000` to reproduce the delivered 528-specification targeted multiverse. The output records the permutation count, nominal P, two-target BH, cross-specification descriptive BH, cell counts and composition for every valid specification.
7. Run script 19 after script 17 to regenerate the fixed-secondary within-sample top-25% KEGG tables and FigSens09. This step filters and visualizes the already computed complete-KEGG result; it does not rerun or alter GSEA statistics.
8. Run script 20 with access to the final Seurat object to reproduce the focused publication-compatibility audit. It retains the paired pseudobulk model and full KEGG BH correction, excludes ad hoc decimal thresholds and q99 caps, and writes the complete result table, target summary and three-patient leading-edge direction audit for the best clean candidate.

Historical scripts contain Windows defaults from the original workstation. The wrappers always pass explicit paths, which should be preferred over those defaults.

## Required R packages

Core packages include Seurat, SeuratObject, Matrix, harmony, edgeR, fgsea, dplyr, tidyr, tibble, ggplot2, patchwork, AnnotationDbi, org.Hs.eg.db, GO.db and infercnv. JAGS is required for the inferCNV installation used here. Exact session records are under environment.

## GitHub availability

The hematologic malignancy project repository contains the code, README, final figures, compact CSV tables and handoff PPT. `03_Raw_Data` and `04_R_Objects` are excluded. Reproduction therefore starts from the frozen compact tables unless the user supplies the GEO raw archive and final R objects locally.
