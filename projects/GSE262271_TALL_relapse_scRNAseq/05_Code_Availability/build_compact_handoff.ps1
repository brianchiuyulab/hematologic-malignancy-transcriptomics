param([Parameter(Mandatory=$true)][string]$HandoffRoot)
$ErrorActionPreference='Stop'
$h=(Resolve-Path -LiteralPath $HandoffRoot).Path
$out=Join-Path $h '00_Quick_Start'
New-Item -ItemType Directory -Path $out -Force | Out-Null
Copy-Item -LiteralPath (Join-Path $h '01_Summary/Leukemia_scRNAseq_STAT5B_reviewed.pptx') -Destination $out -Force
foreach($name in @('Fig01A_cell_lineage_UMAP','Fig01B_annotation_marker_dotplot','Fig05_direction_matched_intersection','Fig06_shared_pathways_dotplot','FigS06_target_definition_comparison')){
 Copy-Item -LiteralPath (Join-Path $h "02_Figures/$name.pdf") -Destination $out -Force
}
foreach($name in @('Table04_scRNAseq_all_KEGG_within15','Table05_bulk_STAT5B_all_KEGG','Table07_same_direction_intersection')){
 Copy-Item -LiteralPath (Join-Path $h "03_Tables/$name.csv") -Destination $out -Force
}
Copy-Item -LiteralPath (Join-Path $h '05_Data_Locations/DATA_AND_OBJECT_LOCATIONS.txt') -Destination $out -Force
# Reader instructions are versioned alongside this script.
Copy-Item -LiteralPath (Join-Path $PSScriptRoot 'QUICK_START.txt') -Destination (Join-Path $out '00_READ_FIRST.txt') -Force
