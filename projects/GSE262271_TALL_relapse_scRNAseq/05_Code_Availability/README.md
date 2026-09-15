# Code availability

Current main analysis: within-sample CNV top15, whole-KEGG FDR < 0.05. R22–26, reviewed 2026-09-15. Deck builder: presentation/build_reviewed_handoff.mjs.

| Script | Purpose |
|---|---|
| R/22_reviewer_audit.R | Marker evidence, cell definitions, count-model refits, high-precision whole-KEGG GSEA, alias sensitivity and paired direction |
| R/23_reviewer_figures.R | Whole-dataset discovery, intersection, target figures and readable tables |
| R/26_shared_pathway_dotplot.R | Aligned40-pathway NES/FDR dot plot and complete20-row presentation pages, called by R23 |
| R/24_independent_reproduction_checks.R | Independent aggregation/model check, historical-ranking comparison, updated lineage Seurat object |
| R/25_reviewed_supporting_figures.R | CNV cutoff distributions and paired gene heatmap |

Use run_reviewed_analysis.ps1 with -ProjectRoot and -BulkRoot. It writes only to an ignored _reproduction_run directory. R23 redraws main figures from reviewed CSV files alone. R25 additionally needs the local reviewed R cache.

R/verify_reviewed_tables.R takes two arguments: the reviewer_audit_results folder
and the 07_Student_Handoff folder. It checks full-collection BH, reproduced
historical ranking, the common tested pathway universe, intersection membership,
directions and sample counts without requiring large R objects.

The paired scRNA-seq model is ~ patient_pair + condition_simple; independent bulk is ~ condition. sc rank = sign(logFC) × sqrt(F); bulk rank = Wald statistic. All tested pathways enter BH. Code checks edgeR4.4.2 and fgsea1.32.4 because another older library exists on the workstation.

See the handoff 04_Code_Availability/REPRODUCIBILITY.txt, reviewer_audit_results, and environment/sessionInfo_reviewer_20260915.txt. Only CSV/text audit files are tracked; large caches are in local 04_R_Objects/reviewed.

## Historical provenance

R01–04 are original preprocessing, annotation and inferCNV code. R05–15 contain earlier KEGG/GO/STAT5/paired/cross-dataset analyses. R16 is the superseded figure renderer. R17–21 retain sensitivity/cutoff analyses. Their output is historical, not the current report. Earlier descriptions of post-hoc settings as pre-specified or publication-defensible must not be reused as prespecification claims.

The original bulk pipeline is included in Bulk_original_pipeline/analysis_gse218858_R.R. It uses the bulk source folder with deposited count files, sample annotation and frozen KEGG2019 Mouse GMT, distinct from the earlier compact deliverables folder.

Historical wrappers support old file layouts. Use the reviewed wrapper for the current outputs. Original frozen results remain in 02_Tables/02_Full_Results for numerical comparison.

## Interpretation

Top15 is the designated main analysis. Because it was fixed after sensitivity analyses, the selection remains exploratory rather than prospectively specified. Pooled reference q95 is not equivalent to within-sample15 or a validated malignancy threshold. Marker review is not independent reference mapping or genotype validation. GSEA on a paired gene ranking is not a paired phenotype-permutation test. Three patient pairs limit generalization. These transcriptomes do not establish drug resistance, mitophagy flux, mitochondrial STAT5 localization or proteasome dependency.

The 15-slide PPT has editable text/tables and embedded R figures. Its builder requires the OpenAI Artifact Tool JavaScript runtime. PNG and vector PDF figures are provided. The complete intersection and Table07 share one display order: increasing max(sc_FDR, bulk_FDR), with alphabetical ties. This maximum is not a combined FDR. Point positions display NES. All pathway names use plain type. Dataset colors are blue (scRNAseq) and orange (bulk). Stars beside the points are computed from unrounded full-catalog BH FDR: *<0.05, **<0.01, ***<0.001. No statistics are recalculated for this display change. Visual layout reference: Reis-de-Oliveira et al., Nature Communications (2024), Figures3B and4D, https://www.nature.com/articles/s41467-024-50875-z .
