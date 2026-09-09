# GSE262271 T-ALL paired diagnosis and relapse scRNA-seq analysis

Clean final package for the paired diagnosis–relapse single-cell RNA-seq analysis of GSE262271.

## Biological question

The analysis asks whether the most CNV-aberrant T-lineage cells at relapse show coordinated transcriptional programs that may connect constitutive STAT5 activity with mitochondrial quality control and proteostasis. The working model is that mitochondrial STAT5 activity could support an oxidative, mitophagy-associated and proteasome-associated state. These data test association and pathway direction. They do not by themselves prove mitochondrial localization, direct STAT5 regulation, drug resistance or therapeutic dependency.

## Primary analysis population

“Potential malignant T cells” are defined as the **CNV-high top 25% within candidate malignant T-lineage cells**, using inferCNV burden estimated against non-T reference populations. This is an operational enrichment strategy, not a definitive clinical malignant-cell call.

- 3 patients with paired diagnosis and relapse samples
- 3,207 cells in the primary CNV-high top-25% population
- Differential expression used patient-level pseudobulk counts
- Paired edgeR design: ~ patient_pair + condition_simple
- GSEA ranking: sign(logFC) × sqrt(F)
- Primary confirmatory threshold: FDR < 0.05
- Exploratory pathway threshold: FDR < 0.10

## Key findings

- Relapse showed strong positive KEGG enrichment for oxidative phosphorylation, ribosome/ribosome biogenesis, DNA replication and spliceosome.
- KEGG mitophagy was relapse-up at NES 1.44, FDR 0.0998. This meets the pre-specified exploratory FDR < 0.10 threshold.
- KEGG proteasome was relapse-up at NES 1.53, FDR 0.109. This is a near-threshold trend, not an FDR < 0.10 hit.
- General autophagy was not significant.
- Custom STAT5 downstream gene sets were directionally positive but not significant.
- In the independent STAT5B N642H bulk RNA-seq dataset, oxidative phosphorylation, mitophagy and proteasome were also positively enriched.
- At KEGG FDR < 0.10, 28 pathways were enriched in the same upward direction in both datasets and none were shared in the downward direction.
- The primary population proportion was heterogeneous across patients and did not show a consistent relapse increase (exact paired Wilcoxon P = 0.75).
- Across 25 predefined sensitivity settings, proteasome and mitophagy NES remained relapse-up in 25/25 settings, but global KEGG FDR < 0.10 was reached in only 10/25 and 6/25 settings, respectively. Direction is robust; threshold-level significance is not.
- Three full-N=3 sensitivity settings placed both pathways below global KEGG FDR 0.10: sample-stratified CNV top 25%, mitochondrial RNA < 4%, and a combined strict-QC rule. These remain sensitivity findings and do not replace the fixed primary analysis.
- A larger targeted multiverse enumerated 528 specifications spanning CNV-high fraction, global versus within-sample CNV thresholds, mitochondrial-RNA and feature cutoffs, upper-feature caps, edgeR filtering/dispersion/normalization and GSEA ranking. Of 520 valid specifications, nominal P < 0.05 occurred in 243 for KEGG mitophagy and 343 for KEGG proteasome; both targets were below 0.05 in 217 specifications. Mitophagy remained relapse-up in 520/520 specifications and proteasome in 518/520. Within-sample CNV thresholds and moderate QC cutoffs generally strengthened the signal, whereas mt < 3% or nFeature >= 1,500 usually weakened it.
- A fixed secondary specification now uses the source-paper QC without extra feature cutoffs and selects the CNV-high top 25% separately within each sample. Across the complete KEGG collection, 62 pathways reached FDR < 0.05 and 87 reached FDR < 0.10. Proteasome (NES 1.60, FDR 0.0779) and mitophagy (NES 1.45, FDR 0.0925) were accompanied by ribosome, oxidative phosphorylation, JAK-STAT, T-cell receptor, cytokine-receptor, ribosome-biogenesis, spliceosome, MAPK, apoptosis and NF-kappa B programs.
- A separate publication-compatibility audit tested 10 pre-specified candidate definitions with full KEGG multiplicity correction while keeping the paired model fixed. The cleanest candidate with both targets below global FDR 0.05 preserved the source-paper QC and selected the CNV-high top 15% within each sample (1,926 cells; minimum 107 cells/sample): proteasome NES 1.89, FDR 0.00586; mitophagy NES 1.69, FDR 0.0152. Leading-edge mean expression increased at relapse in all three patient pairs for both pathways. Adding nFeature >= 500 removed only 26 cells and did not strengthen the result; mt < 4% weakened mitophagy. This remains a transparent exploratory candidate because it was identified through sensitivity analysis and does not retrospectively replace the frozen primary definition.

## Publication figure standard

The final figure suite was rebuilt in September 2026. Diagnosis and relapse use both different colors and different point shapes. Every inferential figure distinguishes FDR < 0.05, exploratory FDR < 0.10 and non-significant results. Gene-level and pathway-level statistics are never mixed. The complete review record is in `01_Figures/FIGURE_AUDIT.csv`.

## Folder map

| Folder | Contents |
|---|---|
| 01_Figures | Main, supplementary and sensitivity figures |
| 02_Tables | Curated key results, complete result tables, sensitivity checks and method inputs |
| 03_Raw_Data | GEO raw archive and source metadata |
| 04_R_Objects | Final Seurat object and inferCNV objects/support files |
| 05_Code_Availability | Exact R scripts, run order, environment records and GitHub notes |
| 06_Handoff_PPT | Lab handoff presentation |

Start with the PPT, then use `01_Figures/FIGURE_INDEX.csv` and `02_Tables/TABLE_GUIDE.csv`. The eight main figures are the recommended reading order. Supplementary and sensitivity figures are separated physically and by filename.

For the large targeted sensitivity grid, start with `multiverse_both_targets_nominal_P_lt_0.05.csv` and `FigSens07_targeted_GSEA_multiverse_specification_curve.png`; the complete one-row-per-pathway-per-specification audit remains available separately.

## Interpretation boundary

The central result is a reproducible pathway-level association between relapse and an oxidative/proteostasis program. The paired sample size is only three patients. Gene-level tests are underpowered and many individual genes are not significant even when GSEA is significant. Functional validation is required to distinguish dependency, adaptation and resistance.

## Data provenance

- Single-cell dataset: GEO GSE262271
- External validation: GSE218858 STAT5B N642H bulk RNA-seq
- Analysis date: July 2026
- Publication-figure and handoff revision: September 2026

This dataset is one analysis module within the broader hematologic malignancy project repository. Large raw and R object files are intentionally excluded from GitHub by `.gitignore`. Their local locations and checksums remain documented in the handoff package.
