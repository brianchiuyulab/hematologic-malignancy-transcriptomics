# Hematologic malignancy transcriptomics

血癌研究的可重現 transcriptomic analyses。請從 [目前交接包](projects/GSE262271_TALL_relapse_scRNAseq/07_Student_Handoff/) 開始。

## GSE262271 and GSE218858

- Human T-ALL scRNA-seq: three diagnosis/relapse patient pairs, paired pseudobulk edgeR.
- Mouse STAT5B N642H bulk RNA-seq: three N642H and three control samples, independent DESeq2.
- Reviewed 2026-09-15: whole-KEGG discovery, direction-matched overlap, target follow-up, conservative lineage annotation and explicit sensitivity definitions.

The designated main analysis uses within-sample CNV top15 and whole-collection KEGG FDR < 0.05. It has 61 up- and 2 down-enriched pathways. Bulk has 107 up and 6 down. Among 275 commonly tested pathway names, 40 are jointly up-enriched and none jointly down-enriched. Alternative cell definitions are supplementary. Top15 was fixed after sensitivity analyses, so its selection remains exploratory.

| Analysis | Proteasome FDR | Mitophagy FDR |
|---|---:|---:|
| Main scRNA-seq within-sample top 15% | 0.008379 | 0.013197 |
| STAT5B N642H bulk | 0.001078 | 0.028234 |
| Sensitivity: original scRNA-seq global top 25% | 0.126157 | 0.126157 |

Both targets have positive NES, not downregulation. Fixed-seed, higher-precision GSEA replaces the old estimates. The original gene statistics reproduce to numerical precision.

Pooled non-T reference q95 and selecting 15% within each sample are different rules. Similar overall fractions do not make top15 unbiased or a validated malignant-cell classifier. Within-run pathway BH does not correct selection among multiple analysis specifications.

## Availability and limits

The repository contains R code, reviewed figures, CSV tables and a 14-slide handoff PPT. Raw GEO archives and large Seurat/inferCNV objects remain local. Their paths and accessions are in the handoff. Original preprocessing and sensitivity results remain available as historical provenance.

These data support transcriptomic association. They do not establish mitochondrial STAT5 localization, mitophagy flux, proteasome dependency or bortezomib resistance.
