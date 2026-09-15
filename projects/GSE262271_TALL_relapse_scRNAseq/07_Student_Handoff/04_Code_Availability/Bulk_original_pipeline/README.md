# GSE218858 R/DESeq2 Code Availability

This folder contains an R-based reproducible analysis for GEO GSE218858.

## Focus

The main comparison is:

`R2S5b_DN` vs `R2b_DN`

This means STAT5B N642H; Rag2-/- DN thymocytes compared with Rag2-/- DN
thymocytes. Other samples are retained only for QC/PCA/context plots.

## Method

- Differential expression: R `DESeq2`
- Ranking: DESeq2 Wald statistic
- GSEA: R `fgsea`
- Gene mapping: `org.Mm.eg.db`
- Gene sets:
  - MSigDB mouse Hallmark
  - MSigDB mouse GO BP
  - Enrichr KEGG_2019_Mouse GMT

## Run

```powershell
Rscript.exe analysis_gse218858_R.R "C:\path\to\R version"
```

The analysis directory must contain:

- `downloads/GSE218858_RAW.tar`
- `downloads/GSE218858_series_matrix.txt.gz`
- `reference/mh.all.v2024.1.Mm.symbols.gmt`
- `reference/m5.go.bp.v2024.1.Mm.symbols.gmt`
- `reference/KEGG_2019_Mouse.enrichr.gmt`

The provided `run_pipeline.ps1` can download these files if they are missing.

