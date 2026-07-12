$ErrorActionPreference = 'Stop'

$repo = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $repo

$env:COMPARE_MAX_TEMPS = '10'
$env:COMPARE_NBOOT = '1'
$env:COMPARE_DROP_RATIO = '0.15'
$env:COMPARE_NSENS = '5'

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
$status = Join-Path $repo 'all_datasets_full_run_status.txt'
"RESUME: $(Get-Date -Format o)" | Add-Content -Path $status

foreach($d in $datasets){
  $base = [System.IO.Path]::GetFileNameWithoutExtension($d).ToLower()
  $out = Join-Path $repo ($base + '_paper_vs_residual_peak_cv_comparison.txt')
  if(Test-Path $out){
    "Skip existing: $d" | Add-Content -Path $status
    continue
  }

  "Running: $d at $(Get-Date -Format o)" | Add-Content -Path $status
  $env:COMPARE_DATASET = $d
  & $matlab -wait -batch run_paper_vs_residual_peak_cv_compare
  if($LASTEXITCODE -ne 0){
    "FAILED: $d exit=$LASTEXITCODE at $(Get-Date -Format o)" | Add-Content -Path $status
    throw "MATLAB failed for $d with exit code $LASTEXITCODE"
  }
  "Done: $d at $(Get-Date -Format o)" | Add-Content -Path $status
}

& .\build_all_datasets_compare_summary.ps1
"Completed: $(Get-Date -Format o)" | Add-Content -Path $status
Write-Host 'Resume run complete.'
