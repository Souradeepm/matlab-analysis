$ErrorActionPreference = 'Stop'

$repo = Split-Path -Parent $MyInvocation.MyCommand.Path
$outCsv = Join-Path $repo 'all_datasets_lambda_values.csv'
$outTxt = Join-Path $repo 'all_datasets_lambda_values_summary.txt'

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

$rows = @()

foreach($d in $datasets){
  $report = Join-Path $repo ($d + '_paper_vs_residual_peak_cv_comparison.txt')
  if(-not (Test-Path $report)){
    continue
  }

  $lines = Get-Content $report
  $headerIdx = -1
  for($i = 0; $i -lt $lines.Count; $i++){
    if($lines[$i] -match '^TemperatureK,ResidualLambda,PaperLambda'){
      $headerIdx = $i
      break
    }
  }

  if($headerIdx -lt 0){
    continue
  }

  for($j = $headerIdx + 1; $j -lt $lines.Count; $j++){
    $line = $lines[$j].Trim()
    if([string]::IsNullOrWhiteSpace($line)){
      continue
    }

    $parts = $line.Split(',')
    if($parts.Count -lt 3){
      continue
    }

    $temp = 0.0; $resLam = 0.0; $paperLam = 0.0
    if(-not [double]::TryParse($parts[0], [ref]$temp)){ continue }
    if(-not [double]::TryParse($parts[1], [ref]$resLam)){ continue }
    if(-not [double]::TryParse($parts[2], [ref]$paperLam)){ continue }

    $rows += [PSCustomObject]@{
      Dataset = $d
      TemperatureK = $temp
      ResidualLambda = $resLam
      PaperLambda = $paperLam
    }
  }
}

$rows | Sort-Object Dataset, TemperatureK | Export-Csv -Path $outCsv -NoTypeInformation -Encoding ascii

$summary = @()
$summary += 'All datasets lambda values summary (15% random-drop CV run)'
$summary += ('Generated: ' + (Get-Date -Format 'yyyy-MM-dd HH:mm:ss'))
$summary += ''
$summary += 'Dataset,Count,ResidualLambdasUnique,PaperLambdasUnique'

foreach($d in $datasets){
  $dsRows = $rows | Where-Object { $_.Dataset -eq $d }
  if($dsRows.Count -eq 0){
    $summary += "$d,0,,"
    continue
  }

  $resUniq = $dsRows.ResidualLambda | Sort-Object -Unique
  $paperUniq = $dsRows.PaperLambda | Sort-Object -Unique

  $resStr = ($resUniq | ForEach-Object { '{0:E4}' -f $_ }) -join '|'
  $paperStr = ($paperUniq | ForEach-Object { '{0:E4}' -f $_ }) -join '|'

  $summary += "$d,$($dsRows.Count),$resStr,$paperStr"
}

Set-Content -Path $outTxt -Value $summary -Encoding ascii
Write-Host "Wrote: $outCsv"
Write-Host "Wrote: $outTxt"
