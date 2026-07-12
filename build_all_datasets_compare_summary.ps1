$ErrorActionPreference = 'Stop'

$repo = Split-Path -Parent $MyInvocation.MyCommand.Path
$out = Join-Path $repo 'all_datasets_paper_vs_residual_peak_cv_comparison.txt'

$datasets = @(
  's2022sap',
  's2022al',
  's2222sap',
  's2222al',
  's2302sap',
  's2302al',
  's2322sap',
  's2322al',
  's2332sap',
  's2422sap',
  's2422al'
)

$lines = @()
$lines += 'All-dataset paper vs residual comparison (15% drop residual-change sensitivity)'
$lines += ('Generated: ' + (Get-Date -Format 'yyyy-MM-dd HH:mm:ss'))
$lines += ''
$lines += 'Dataset,Status,Temps,MeanResidualPeaks,MeanPaperPeaks,MeanPeakDelta,MeanResidualSensitivity,MeanPaperSensitivity,MeanSensitivityDelta,PaperLowerSensitivityCount,PaperHigherSensitivityCount,ReportFile'

foreach($d in $datasets){
  $f = Join-Path $repo ($d + '_paper_vs_residual_peak_cv_comparison.txt')
  if(-not (Test-Path $f)){
    $lines += "$d,missing,,,,,,,,,,$f"
    continue
  }

  $txt = Get-Content $f

  function Get-Val([string]$pattern){
    $m = $txt | Select-String -Pattern $pattern | Select-Object -First 1
    if($null -eq $m){ return '' }
    return ($m.Line -replace '^[^:]+:\s*','').Trim()
  }

  $temps = Get-Val '^Total temperatures:'
  $mrp = Get-Val '^Mean peaks \(residual method\):'
  $mpp = Get-Val '^Mean peaks \(paper method\)\s*:'
  $mrs = Get-Val '^Mean residual-change sensitivity % \(residual-selected lambda\):'
  $mps = Get-Val '^Mean residual-change sensitivity % \(paper-selected lambda\)\s*:'
  $plc = Get-Val '^Temps paper lower sensitivity\s*:'
  $phc = Get-Val '^Temps paper higher sensitivity\s*:'

  $mpd = ''
  $msd = ''
  $mrpNum = 0.0; $mppNum = 0.0; $mrsNum = 0.0; $mpsNum = 0.0
  if([double]::TryParse($mrp, [ref]$mrpNum) -and [double]::TryParse($mpp, [ref]$mppNum)){
    $mpd = ('{0:F6}' -f [Math]::Abs($mppNum - $mrpNum))
  }
  if([double]::TryParse($mrs, [ref]$mrsNum) -and [double]::TryParse($mps, [ref]$mpsNum)){
    $msd = ('{0:E8}' -f [Math]::Abs($mpsNum - $mrsNum))
  }

  $lines += "$d,ok,$temps,$mrp,$mpp,$mpd,$mrs,$mps,$msd,$plc,$phc,$f"
}

Set-Content -Path $out -Value $lines -Encoding ascii
Write-Host "Wrote: $out"
