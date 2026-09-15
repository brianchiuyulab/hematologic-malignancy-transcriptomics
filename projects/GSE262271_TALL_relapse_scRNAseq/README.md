# T-ALL relapse and STAT5B transcriptomics

Current handoff: [07_Student_Handoff/00_READ_ME.txt](07_Student_Handoff/00_READ_ME.txt).

The 14-slide PPT follows datasets, lineage annotation, paired model, whole-KEGG enrichment in each dataset, same-direction intersection, target pathways, limitations and reproduction.

## Reviewed results, 2026-09-15

| Dataset / definition | Up, FDR < 0.10 | Down, FDR < 0.10 |
|---|---:|---:|
| scRNA-seq, exploratory within-sample top15 | 83 | 5 |
| scRNA-seq, original global top25 | 47 | 1 |
| STAT5B N642H bulk RNA-seq | 128 | 8 |

Among 275 commonly tested pathway names, top15 and bulk share 54 up-enriched pathways and zero down-enriched pathways. Original global25 and bulk share 23 up-enriched pathways.

Proteasome and Mitophagy are up-enriched in top15 (FDR 0.008379 and 0.013197) and bulk (0.001078 and 0.028234). Neither passes FDR < 0.10 in original global25 after high-precision GSEA (both 0.126157). General autophagy is not significant.

## Corrections and limits

The old cross-dataset figure used global25 while later significant results used within-sample15. Definitions are now visible. Gene and pathway significance are separate. Original gene statistics reproduce using edgeR4.4.2 and the stored count object; bulk DESeq2 statistics also reproduce. Fixed-seed high-precision GSEA replaces unstable near-threshold estimates.

Top15 remains exploratory. A pooled reference q95=0.0269455, exceeded by 15.419% of T-lineage candidates, does not validate taking exactly15% within each sample. Non-T reference cells are not independently verified healthy-cell ground truth. Pathway BH does not correct the search over specifications.

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
| 05_Code_Availability | R22–25, original pipeline, environment and audit records |

Current renderer: R23, supporting figures: R25. Legacy R16 and old decks are not current. Use run_reviewed_analysis.ps1 for isolated reproduction. This audit did not realign FASTQ or rerun Cell Ranger, Harmony or inferCNV.

The findings do not prove mitochondrial STAT5 localization, mitophagy activity, proteasome dependency or bortezomib resistance.
