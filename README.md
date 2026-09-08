# Hematologic malignancy transcriptomics

這個 repository 用來集中血液腫瘤的可重現 transcriptomic analyses。每個 dataset 都放在 `projects/` 下，使用一致的交接結構：final figures、compact tables、code availability、handoff PowerPoint，以及只保留說明檔的 raw-data / R-object folders。

## Current project modules

### GSE262271: paired T-ALL diagnosis and relapse scRNA-seq

- Primary population: CNV-high top 25% within candidate malignant T-lineage cells, operationally labeled **potential malignant T cells**.
- Statistical design: patient-level paired pseudobulk edgeR, `~ patient_pair + condition_simple`.
- Main result: relapse-associated cells show coordinated ribosome / ribosome-biogenesis, spliceosome, DNA-replication and oxidative-phosphorylation programs.
- Hypothesis-focused result: KEGG mitophagy is exploratory at FDR 0.0998; proteasome is near threshold at FDR 0.109 and is not counted as FDR < 0.10.
- STAT5 boundary: STAT5A is nominally increased, but gene-level FDR and tested STAT5 downstream gene sets are not significant.
- External comparison: STAT5B N642H bulk RNA-seq shares upward OXPHOS, mitophagy, proteasome and ribosome directions, supporting association rather than direct causality.

Start here: [`projects/GSE262271_TALL_relapse_scRNAseq/START_HERE_請先看.txt`](projects/GSE262271_TALL_relapse_scRNAseq/START_HERE_%E8%AB%8B%E5%85%88%E7%9C%8B.txt)

## Repository policy

- GitHub contains publication figures, compact CSV tables, analysis code, environment records and the reviewed handoff deck.
- GEO archives, large Seurat / inferCNV objects and reproduction workspaces remain local and are ignored.
- Main figures use colorblind-aware palettes, explicit direction / evidence symbols and paired-patient labeling.
- FDR < 0.05 is the primary significance threshold. FDR < 0.10 is labeled exploratory. Near-threshold results are reported numerically and are not promoted to significant.

## Directory layout

```text
projects/
  GSE262271_TALL_relapse_scRNAseq/
    01_Figures/
    02_Tables/
    03_Raw_Data/          # documentation only on GitHub
    04_R_Objects/         # documentation only on GitHub
    05_Code_Availability/
    06_Handoff_PPT/
```

## Interpretation boundary

The analyses establish paired transcriptomic associations and pathway direction. They do not by themselves prove mitochondrial STAT5 localization, mitophagy flux, proteasome dependency or bortezomib resistance. Those claims require functional validation.
