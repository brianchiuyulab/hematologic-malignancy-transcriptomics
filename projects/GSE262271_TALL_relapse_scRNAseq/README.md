# T-ALL relapse and STAT5B transcriptomics

Current handoff: [07_Student_Handoff/00_READ_ME.txt](07_Student_Handoff/00_READ_ME.txt).

The 15-slide PPT follows datasets, lineage annotation, paired model, whole-KEGG enrichment in each dataset, same-direction intersection, target pathways, limitations and reproduction. Figure06 shows all40 shared pathways as aligned NES dot plots; slides9–10 show the same order in two readable pages. No pathways are omitted from the intersection presentation.

## Reviewed results, 2026-09-15

| Dataset / definition | Up, FDR < 0.05 | Down, FDR < 0.05 |
|---|---:|---:|
| Main scRNA-seq, within-sample top15 | 61 | 2 |
| Sensitivity scRNA-seq, original global top25 | 28 | 0 |
| STAT5B N642H bulk RNA-seq | 107 | 6 |

Among 275 commonly tested pathway names, top15 and bulk share 40 up-enriched pathways and zero down-enriched pathways. Original global25 and bulk share 13 up-enriched pathways.

Proteasome and Mitophagy are up-enriched in top15 (FDR 0.008379 and 0.013197) and bulk (0.001078 and 0.028234). Neither passes FDR < 0.05 in original global25 after high-precision GSEA (both 0.126157). General autophagy is not significant.

## Corrections and limits

The old cross-dataset figure used global25 while later significant results used within-sample15. Definitions are now visible. Gene and pathway significance are separate. Original gene statistics reproduce using edgeR4.4.2 and the stored count object; bulk DESeq2 statistics also reproduce. Fixed-seed high-precision GSEA replaces unstable near-threshold estimates.

Within-sample top15 is the designated main analysis and FDR < 0.05 is the common reporting threshold. Its selection after sensitivity analysis remains exploratory. A pooled reference q95=0.0269455, exceeded by 15.419% of T-lineage candidates, does not uniquely determine taking exactly15% within each sample. Only 1,407 of the 1,926 top15 cells (73.1%) also exceed the fixed q95 score; sample-level comparisons are in TableS03 and TableS06. Non-T reference cells are not independently verified healthy-cell ground truth. Pathway BH does not correct the search over specifications.

Broad lineage labels have cluster-marker support. HSC-like/quiescent/ER-stress are no longer presented as validated T-ALL subtypes. CNV-high potential malignant status is a separate operational selection. Independent malignancy and same-sample doublet validation remain needed.

## Folder map

| Folder | Purpose |
|---|---|
| 07_Student_Handoff | One reader package: PPT, figures, tables, code, data locations |
| 01_Figures | Current reviewed figures only |
| 02_Tables/01_Key_Results | Current readable CSV tables |
| 02_Tables/02_Full_Results | Historical frozen model/GSEA inputs for reproduction |
| 02_Tables/03_Sensitivity_Checks | Historical sensitivity results, separate from current results |
| 02_Tables/04_Method_Inputs | Original cell metadata / CNV inputs |
| 03_Raw_Data | GEO archive, not uploaded |
| 04_R_Objects | Original and reviewed Seurat / inferCNV objects, not uploaded |
| 05_Code_Availability | R22–26, original pipeline, environment and audit records |

Current renderer: R23, supporting figures: R25. Legacy R16 and old decks are not current. Use run_reviewed_analysis.ps1 for isolated reproduction. This audit did not realign FASTQ or rerun Cell Ranger, Harmony or inferCNV.

The findings do not prove mitochondrial STAT5 localization, mitophagy activity, proteasome dependency or bortezomib resistance.
