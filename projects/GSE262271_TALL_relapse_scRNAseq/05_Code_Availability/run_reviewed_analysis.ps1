param(
  [Parameter(Mandatory=$true)][string]$ProjectRoot,
  [Parameter(Mandatory=$true)][string]$BulkRoot,
  [string]$Rscript = 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe',
  [string]$RLibrary = 'C:\Users\User\Documents\R\win-library\4.4'
)
$ErrorActionPreference = 'Stop'
$ProjectRoot = (Resolve-Path -LiteralPath $ProjectRoot).Path
$BulkRoot = (Resolve-Path -LiteralPath $BulkRoot).Path
$run = Join-Path $ProjectRoot '_reproduction_run\reviewed_20260915'
$results = Join-Path $run 'results'
$delivery = Join-Path $run 'delivery'
New-Item -ItemType Directory -Path $results -Force | Out-Null
$previousLocale = $env:LC_ALL
$previousLibrary = $env:GSE_REPRO_R_LIB
$previousResume = $env:AUDIT_RESUME
$previousGseaReuse = $env:AUDIT_GSEA_REUSE
try {
  $env:LC_ALL = 'English_United States.utf8'
  $env:GSE_REPRO_R_LIB = $RLibrary
  $env:AUDIT_RESUME = '0'
  $env:AUDIT_GSEA_REUSE = '0'
  & $Rscript (Join-Path $PSScriptRoot 'R\22_reviewer_audit.R') $ProjectRoot $BulkRoot $results
  if ($LASTEXITCODE -ne 0) { throw 'Count-model/GSEA audit failed.' }
  & $Rscript (Join-Path $PSScriptRoot 'R\23_reviewer_figures.R') $results $delivery
  if ($LASTEXITCODE -ne 0) { throw 'Main figures failed.' }
  & $Rscript (Join-Path $PSScriptRoot 'R\24_independent_reproduction_checks.R') $ProjectRoot $results
  if ($LASTEXITCODE -ne 0) { throw 'Independent reproduction checks failed.' }
  & $Rscript (Join-Path $PSScriptRoot 'R\25_reviewed_supporting_figures.R') $results $delivery
  if ($LASTEXITCODE -ne 0) { throw 'Supporting figures failed.' }
  Write-Output "Reproduction complete: $delivery"
} finally {
  $env:LC_ALL = $previousLocale
  $env:GSE_REPRO_R_LIB = $previousLibrary
  $env:AUDIT_RESUME = $previousResume
  $env:AUDIT_GSEA_REUSE = $previousGseaReuse
}
