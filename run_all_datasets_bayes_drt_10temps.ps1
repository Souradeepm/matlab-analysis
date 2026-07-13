$ErrorActionPreference = 'Stop'

$repo = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $repo

$env:BAYES_DRT_MAX_TEMPS = '10'
$env:BAYES_DRT_LAMBDA_MIN_EXP = '-8'
$env:BAYES_DRT_LAMBDA_MAX_EXP = '2'
$env:BAYES_DRT_LAMBDA_COUNT = '21'

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
$status = Join-Path $repo 'all_datasets_bayes_drt_10temps_status.txt'
$summary = Join-Path $repo 'all_datasets_bayes_drt_10temps_summary.txt'
"Start: $(Get-Date -Format o)" | Set-Content -Path $status -Encoding ascii

$summaryRows = @()

Get-ChildItem -Path $repo -Filter '*_bayes_drt_matlab2011_10temp.txt' -File -ErrorAction SilentlyContinue |
  Remove-Item -Force -ErrorAction SilentlyContinue

foreach($d in $datasets){
  ('Running dataset: {0} at {1}' -f $d, (Get-Date -Format o)) | Add-Content -Path $status
  $env:BAYES_DRT_DATASET = $d

  $ok = $false
  for($attempt = 1; $attempt -le 2; $attempt++){
    '  Attempt {0}' -f $attempt | Add-Content -Path $status

    Get-Process -Name MATLAB -ErrorAction SilentlyContinue | Stop-Process -Force -ErrorAction SilentlyContinue

    $outLog = Join-Path $repo ([System.IO.Path]::GetFileNameWithoutExtension($d).ToLower() + '_bayes_drt_10temp_attempt' + $attempt + '_out.log')
    $errLog = Join-Path $repo ([System.IO.Path]::GetFileNameWithoutExtension($d).ToLower() + '_bayes_drt_10temp_attempt' + $attempt + '_err.log')
    $reportPath = Join-Path $repo ([System.IO.Path]::GetFileNameWithoutExtension($d).ToLower() + '_bayes_drt_matlab2011_10temp.txt')
    $p = Start-Process -FilePath $matlab -ArgumentList '-wait','-batch','run_bayes_drt_workflow_matlab2011' -PassThru -NoNewWindow -RedirectStandardOutput $outLog -RedirectStandardError $errLog
    $p.WaitForExit()
    $code = $p.ExitCode
    if($null -eq $code -or $code -eq ''){
      if(Test-Path $reportPath){
        $code = 0
      } else {
        $code = -998
      }
    }
    if($code -eq 0 -and (Test-Path $reportPath)){
      $ok = $true
      "  Attempt $attempt success" | Add-Content -Path $status
      break
    }

    "  Attempt $attempt failed exit=$code" | Add-Content -Path $status
  }

  if(-not $ok){
    "FAILED: $d after retries at $(Get-Date -Format o)" | Add-Content -Path $status
    continue
  }

  $reportPath = Join-Path $repo ([System.IO.Path]::GetFileNameWithoutExtension($d).ToLower() + '_bayes_drt_matlab2011_10temp.txt')
  if(-not (Test-Path $reportPath)){
    "FAILED to parse report for $d" | Add-Content -Path $status
    continue
  }

  $lines = Get-Content $reportPath
  $parsed = [ordered]@{}
  foreach($line in $lines){
    if($line -match '^Mean selected lambda:\s+([0-9.eE+-]+)') { $parsed.MeanSelectedLambda = [double]$Matches[1]; continue }
    if($line -match '^Mean real-part CV:\s+([0-9.eE+-]+)') { $parsed.MeanRealCV = [double]$Matches[1]; continue }
    if($line -match '^Mean imag-part CV:\s+([0-9.eE+-]+)') { $parsed.MeanImagCV = [double]$Matches[1]; continue }
    if($line -match '^Mean total CV:\s+([0-9.eE+-]+)') { $parsed.MeanTotalCV = [double]$Matches[1]; continue }
    if($line -match '^Mean absolute residual:\s+([0-9.eE+-]+)') { $parsed.MeanResidual = [double]$Matches[1]; continue }
    if($line -match '^Mean peak count:\s+([0-9.eE+-]+)') { $parsed.MeanPeakCount = [double]$Matches[1]; continue }
    if($line -match '^Total processed temperatures:\s+([0-9]+)') { $parsed.TemperatureCount = [int]$Matches[1]; continue }
  }
  $parsed = [PSCustomObject]$parsed

  $summaryRows += [PSCustomObject]@{
    Dataset = [System.IO.Path]::GetFileNameWithoutExtension($d).ToLower()
    Temperatures = $parsed.TemperatureCount
    MeanSelectedLambda = $parsed.MeanSelectedLambda
    MeanRealCV = $parsed.MeanRealCV
    MeanImagCV = $parsed.MeanImagCV
    MeanTotalCV = $parsed.MeanTotalCV
    MeanResidual = $parsed.MeanResidual
    MeanPeakCount = $parsed.MeanPeakCount
  }

  "Done dataset: $d at $(Get-Date -Format o)" | Add-Content -Path $status
}

$summaryLines = @()
$summaryLines += 'Dataset,Temperatures,MeanSelectedLambda,MeanRealCV,MeanImagCV,MeanTotalCV,MeanResidual,MeanPeakCount'
foreach($row in $summaryRows){
  $summaryLines += ('{0},{1},{2:E8},{3:E8},{4:E8},{5:E8},{6:E8},{7:F4}' -f `
    $row.Dataset, $row.Temperatures, $row.MeanSelectedLambda, $row.MeanRealCV, $row.MeanImagCV, $row.MeanTotalCV, $row.MeanResidual, $row.MeanPeakCount)
}
$summaryLines | Set-Content -Path $summary -Encoding ascii
"Completed: $(Get-Date -Format o)" | Add-Content -Path $status
Write-Host 'Bayes-DRT 10-temperature workflow complete.'
