param(
    [Parameter(Mandatory=$true)]
    [string]$AnalysisDir,

    [string]$Rscript = "C:\Program Files\R\R-4.4.2\bin\Rscript.exe"
)

$ErrorActionPreference = "Stop"

$analysisDirFull = [IO.Path]::GetFullPath($AnalysisDir)
$downloadDir = Join-Path $analysisDirFull "downloads"
$referenceDir = Join-Path $analysisDirFull "reference"

New-Item -ItemType Directory -Force -Path $downloadDir | Out-Null
New-Item -ItemType Directory -Force -Path $referenceDir | Out-Null

function Download-IfMissing {
    param(
        [string]$Url,
        [string]$Destination
    )
    if ((-not (Test-Path -LiteralPath $Destination)) -or ((Get-Item -LiteralPath $Destination).Length -eq 0)) {
        Write-Host "Downloading $Url"
        Invoke-WebRequest -Uri $Url -OutFile $Destination -UseBasicParsing
    } else {
        Write-Host "Already exists $Destination"
    }
}

$geoBase = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE218nnn/GSE218858"
Download-IfMissing "$geoBase/suppl/GSE218858_RAW.tar" (Join-Path $downloadDir "GSE218858_RAW.tar")
Download-IfMissing "$geoBase/matrix/GSE218858_series_matrix.txt.gz" (Join-Path $downloadDir "GSE218858_series_matrix.txt.gz")
Download-IfMissing "$geoBase/suppl/GSE218858_all_DIFF_N642H_DE_result_tables.xlsx" (Join-Path $downloadDir "GSE218858_all_DIFF_N642H_DE_result_tables.xlsx")
Download-IfMissing "$geoBase/suppl/GSE218858_all_GT_N642H_DE_result_tables.xlsx" (Join-Path $downloadDir "GSE218858_all_GT_N642H_DE_result_tables.xlsx")

Download-IfMissing "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Mm/mh.all.v2024.1.Mm.symbols.gmt" (Join-Path $referenceDir "mh.all.v2024.1.Mm.symbols.gmt")
Download-IfMissing "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Mm/m5.go.bp.v2024.1.Mm.symbols.gmt" (Join-Path $referenceDir "m5.go.bp.v2024.1.Mm.symbols.gmt")
Download-IfMissing "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Mouse" (Join-Path $referenceDir "KEGG_2019_Mouse.enrichr.gmt")

& $Rscript (Join-Path $PSScriptRoot "analysis_gse218858_R.R") $analysisDirFull
