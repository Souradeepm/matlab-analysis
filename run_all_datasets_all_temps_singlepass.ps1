$ErrorActionPreference = 'Continue'

$repo = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $repo

$env:COMPARE_MAX_TEMPS = '100000'
$env:COMPARE_NBOOT = '1'
$env:COMPARE_DROP_RATIO = '0.15'
$env:COMPARE_NSENS = '1'

$datasets = @(
  'S2022Sap.xlsx',
  'S2022Al.xlsx',
  'S2222Sap.xlsx',
  'S2222Al.xlsx',
  'S2302Sap.xlsx',
  'S2302Al.xlsx',
  'S2322Sap.xlsx',
  'S2322Al.xlsx',
  'S2332Sap.xlsx',
  'S2422Sap.xlsx',
  'S2422Al.xlsx'
)

$matlab = 'C:\Progra~1\MATLAB\R2023a\bin\matlab.exe'
$status = Join-Path $repo 'all_datasets_all_temps_status_clean.txt'
"Start: $(Get-Date -Format o)" | Set-Content -Path $status -Encoding ascii

Get-ChildItem -Path $repo -Filter '*_paper_vs_residual_peak_cv_comparison.txt' -File -ErrorAction SilentlyContinue |
  Where-Object { $_.Name -ne 'all_datasets_paper_vs_residual_peak_cv_comparison.txt' } |
  Remove-Item -Force -ErrorAction SilentlyContinue

foreach($d in $datasets){
  "Running: $d at $(Get-Date -Format o)" | Add-Content -Path $status
  $env:COMPARE_DATASET = $d

  Get-Process -Name MATLAB -ErrorAction SilentlyContinue | Stop-Process -Force -ErrorAction SilentlyContinue

  & $matlab -wait -batch run_paper_vs_residual_peak_cv_compare
  $code = $LASTEXITCODE

  if($code -eq 0){
    "Done: $d at $(Get-Date -Format o)" | Add-Content -Path $status
  } else {
    "FAILED: $d exit=$code at $(Get-Date -Format o)" | Add-Content -Path $status
  }
}

& .\build_all_datasets_compare_summary.ps1
& .\build_all_datasets_lambda_values.ps1
"Completed: $(Get-Date -Format o)" | Add-Content -Path $status
Write-Host 'Single-pass all-datasets all-temperatures run finished.'
