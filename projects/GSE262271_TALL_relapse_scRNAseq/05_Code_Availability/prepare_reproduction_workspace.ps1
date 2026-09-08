param(
  [string]$ProjectRoot = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
)

$ErrorActionPreference = 'Stop'
$work = Join-Path $ProjectRoot '_reproduction_run'
$output = Join-Path $work 'outputs/GSE262271_TALL_relapse_scRNAseq_final'

$dirs = @(
  'raw_data', 'source_metadata', 'processed_objects', 'infercnv_inputs',
  'results', 'figures', 'logs',
  'outputs/GSE262271_TALL_relapse_scRNAseq_final/csv/potential_malignant_T_cells',
  'outputs/GSE262271_TALL_relapse_scRNAseq_final/csv/sensitivity_checks',
  'outputs/GSE262271_TALL_relapse_scRNAseq_final/png/potential_malignant_T_cells',
  'outputs/GSE262271_TALL_relapse_scRNAseq_final/png/sensitivity_checks'
)
foreach ($dir in $dirs) {
  New-Item -ItemType Directory -Path (Join-Path $work $dir) -Force | Out-Null
}

Copy-Item -LiteralPath (Join-Path $ProjectRoot '03_Raw_Data/GSE262271_RAW.tar') -Destination (Join-Path $work 'raw_data') -Force
Copy-Item -Path (Join-Path $ProjectRoot '03_Raw_Data/source_metadata/*') -Destination (Join-Path $work 'source_metadata') -Force
Copy-Item -LiteralPath (Join-Path $ProjectRoot '04_R_Objects/GSE262271_seurat_qc_harmony_annotated.rds') -Destination (Join-Path $work 'processed_objects') -Force
Copy-Item -Path (Join-Path $ProjectRoot '04_R_Objects/inferCNV/*') -Destination (Join-Path $work 'infercnv_inputs') -Force

$csvTarget = Join-Path $output 'csv/potential_malignant_T_cells'
Copy-Item -Path (Join-Path $ProjectRoot '02_Tables/01_Key_Results/*.csv') -Destination $csvTarget -Force
Copy-Item -Path (Join-Path $ProjectRoot '02_Tables/02_Full_Results/*.csv') -Destination $csvTarget -Force
Copy-Item -Path (Join-Path $ProjectRoot '02_Tables/04_Method_Inputs/*.csv') -Destination $csvTarget -Force
Copy-Item -Path (Join-Path $ProjectRoot '02_Tables/03_Sensitivity_Checks/*.csv') -Destination (Join-Path $output 'csv/sensitivity_checks') -Force

Write-Host "Prepared isolated reproduction workspace: $work"
