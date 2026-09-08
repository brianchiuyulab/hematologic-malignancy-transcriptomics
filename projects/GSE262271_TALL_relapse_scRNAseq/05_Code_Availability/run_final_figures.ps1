param(
  [string]$ProjectRoot = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path,
  [string]$Rscript = 'Rscript'
)

$ErrorActionPreference = 'Stop'
$output = Join-Path $ProjectRoot '_reproduction_run/publication_figures'
$r = Join-Path $PSScriptRoot 'R'

& $Rscript (Join-Path $r '16_publication_figure_suite.R') $ProjectRoot $output
if ($LASTEXITCODE -ne 0) { throw 'Publication figure rendering failed.' }

Write-Host "Recreated publication figures under: $output"
