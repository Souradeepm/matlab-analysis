$ErrorActionPreference = 'Stop'

$repo = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $repo

$env:COMPARE_NBOOT = '1'
$env:COMPARE_DROP_RATIO = '0.15'
$env:COMPARE_NSENS = '1'
$env:COMPARE_MAX_TEMPS = '10'

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
$status = Join-Path $repo 'all_datasets_all_temps_batched10_status.txt'
"Start: $(Get-Date -Format o)" | Set-Content -Path $status -Encoding ascii
$logDir = Join-Path $repo 'batch_logs'
$batchTimeoutSec = 420
$tmp = 0.0
if([double]::TryParse($env:COMPARE_BATCH_TIMEOUT_SEC, [ref]$tmp) -and $tmp -ge 60){
  $batchTimeoutSec = [int][Math]::Round($tmp)
}
$cleanOutputs = $false
if([double]::TryParse($env:COMPARE_CLEAN_OUTPUTS, [ref]$tmp) -and $tmp -ge 1){
  $cleanOutputs = $true
}
if(-not (Test-Path $logDir)){
  New-Item -ItemType Directory -Path $logDir | Out-Null
}

# Clean previous logs and outputs only when explicitly requested.
if($cleanOutputs){
  Get-ChildItem -Path $logDir -Filter '*.log' -File -ErrorAction SilentlyContinue |
    Remove-Item -Force -ErrorAction SilentlyContinue

  Get-ChildItem -Path $repo -Filter '*_paper_vs_residual_peak_cv_comparison_b*.txt' -File -ErrorAction SilentlyContinue |
    Remove-Item -Force -ErrorAction SilentlyContinue
}

function Merge-DatasetBatchReports {
  param([string]$DatasetFile)

  $base = [System.IO.Path]::GetFileNameWithoutExtension($DatasetFile).ToLower()
  $pattern = "$base`_paper_vs_residual_peak_cv_comparison_b*.txt"
  $files = Get-ChildItem -Path $repo -Filter $pattern -File -ErrorAction SilentlyContinue |
    Sort-Object {
      $m = [regex]::Match($_.BaseName, 'b(\d+)_')
      if($m.Success){ [int]$m.Groups[1].Value } else { 999999 }
    }

  if($files.Count -eq 0){
    return
  }

  $rows = @()
  foreach($f in $files){
    $lines = Get-Content $f.FullName
    $headerIdx = -1
    for($i=0; $i -lt $lines.Count; $i++){
      if($lines[$i] -match '^TemperatureK,ResidualLambda,PaperLambda'){
        $headerIdx = $i
        break
      }
    }
    if($headerIdx -lt 0){ continue }

    for($j=$headerIdx+1; $j -lt $lines.Count; $j++){
      $line = $lines[$j].Trim()
      if([string]::IsNullOrWhiteSpace($line)){ continue }
      $parts = $line.Split(',')
      if($parts.Count -lt 9){ continue }

      $t=0.0; $rl=0.0; $pl=0.0; $rp=0; $pp=0; $rs=0.0; $ps=0.0
      if(-not [double]::TryParse($parts[0], [ref]$t)){ continue }
      if(-not [double]::TryParse($parts[1], [ref]$rl)){ continue }
      if(-not [double]::TryParse($parts[2], [ref]$pl)){ continue }
      if(-not [int]::TryParse($parts[3], [ref]$rp)){ continue }
      if(-not [int]::TryParse($parts[4], [ref]$pp)){ continue }
      if(-not [double]::TryParse($parts[6], [ref]$rs)){ continue }
      if(-not [double]::TryParse($parts[7], [ref]$ps)){ continue }

      $rows += [PSCustomObject]@{
        TemperatureK = $t
        ResidualLambda = $rl
        PaperLambda = $pl
        ResidualPeaks = $rp
        PaperPeaks = $pp
        ResidualSensitivity = $rs
        PaperSensitivity = $ps
      }
    }
  }

  if($rows.Count -eq 0){
    return
  }

  $rows = $rows | Sort-Object TemperatureK

  $deltaPeakSigned = @(); $deltaPeakAbs = @(); $deltaSensSigned = @(); $deltaSensAbs = @()
  foreach($r in $rows){
    $dp = [double]$r.PaperPeaks - [double]$r.ResidualPeaks
    $ds = [double]$r.PaperSensitivity - [double]$r.ResidualSensitivity
    $deltaPeakSigned += $dp
    $deltaPeakAbs += [Math]::Abs($dp)
    $deltaSensSigned += $ds
    $deltaSensAbs += [Math]::Abs($ds)
  }

  $meanResPeaks = ($rows | Measure-Object -Property ResidualPeaks -Average).Average
  $meanPaperPeaks = ($rows | Measure-Object -Property PaperPeaks -Average).Average
  $meanPeakDeltaAbs = ($deltaPeakAbs | Measure-Object -Average).Average
  $medianPeakDeltaAbs = ($deltaPeakAbs | Sort-Object)[[int][Math]::Floor(($deltaPeakAbs.Count-1)/2)]

  $meanResSens = ($rows | Measure-Object -Property ResidualSensitivity -Average).Average
  $meanPaperSens = ($rows | Measure-Object -Property PaperSensitivity -Average).Average
  $meanSensDeltaAbs = ($deltaSensAbs | Measure-Object -Average).Average
  $medianSensDeltaAbs = ($deltaSensAbs | Sort-Object)[[int][Math]::Floor(($deltaSensAbs.Count-1)/2)]

  $morePeaks = ($deltaPeakSigned | Where-Object { $_ -gt 0 }).Count
  $equalPeaks = ($deltaPeakSigned | Where-Object { $_ -eq 0 }).Count
  $fewerPeaks = ($deltaPeakSigned | Where-Object { $_ -lt 0 }).Count

  $lowerSens = ($deltaSensSigned | Where-Object { $_ -lt 0 }).Count
  $equalSens = ($deltaSensSigned | Where-Object { $_ -eq 0 }).Count
  $higherSens = ($deltaSensSigned | Where-Object { $_ -gt 0 }).Count

  $out = Join-Path $repo ($base + '_paper_vs_residual_peak_cv_comparison.txt')
  $o = @()
  $o += 'Paper vs minimum-residual DRT comparison'
  $o += ('Dataset: ' + $DatasetFile)
  $o += ('Generated: ' + (Get-Date -Format 'dd-MMM-yyyy HH:mm:ss'))
  $o += ('Total temperatures: ' + $rows.Count)
  $o += ''
  $o += 'Bootstrap repeats for variance metric: 1'
  $o += ''
  $o += 'Drop ratio for sensitivity metric: 0.15'
  $o += 'Sensitivity repeats per lambda: 1'
  $o += ''
  $o += 'Aggregate comparison'
  $o += ('Mean peaks (residual method): {0:F4}' -f $meanResPeaks)
  $o += ('Mean peaks (paper method)   : {0:F4}' -f $meanPaperPeaks)
  $o += ('Mean peak delta abs |paper-res| : {0:F4}' -f $meanPeakDeltaAbs)
  $o += ('Median peak delta abs |paper-res|: {0:F4}' -f $medianPeakDeltaAbs)
  $o += ('Temps paper has more peaks      : {0}' -f $morePeaks)
  $o += ('Temps equal peaks               : {0}' -f $equalPeaks)
  $o += ('Temps paper has fewer peaks     : {0}' -f $fewerPeaks)
  $o += ''
  $o += ('Mean residual-change sensitivity % (residual-selected lambda): {0:E8}' -f $meanResSens)
  $o += ('Mean residual-change sensitivity % (paper-selected lambda)   : {0:E8}' -f $meanPaperSens)
  $o += ('Mean sensitivity delta abs % |paper-res|                    : {0:E8}' -f $meanSensDeltaAbs)
  $o += ('Median sensitivity delta abs % |paper-res|                  : {0:E8}' -f $medianSensDeltaAbs)
  $o += ('Temps paper lower sensitivity                             : {0}' -f $lowerSens)
  $o += ('Temps equal sensitivity                                   : {0}' -f $equalSens)
  $o += ('Temps paper higher sensitivity                            : {0}' -f $higherSens)
  $o += ''
  $o += 'Per-temperature detail'
  $o += 'TemperatureK,ResidualLambda,PaperLambda,ResidualPeaks,PaperPeaks,PeakDeltaAbs_PaperMinusRes,ResidualSensitivity,PaperSensitivity,SensitivityDeltaAbs_PaperMinusRes'

  foreach($r in $rows){
    $dp = [Math]::Abs([double]$r.PaperPeaks - [double]$r.ResidualPeaks)
    $ds = [Math]::Abs([double]$r.PaperSensitivity - [double]$r.ResidualSensitivity)
    $o += ('{0:F6},{1:E8},{2:E8},{3},{4},{5},{6:E8},{7:E8},{8:E8}' -f $r.TemperatureK,$r.ResidualLambda,$r.PaperLambda,$r.ResidualPeaks,$r.PaperPeaks,$dp,$r.ResidualSensitivity,$r.PaperSensitivity,$ds)
  }

  Set-Content -Path $out -Value $o -Encoding ascii
}

foreach($d in $datasets){
  ('Running dataset: {0} at {1}' -f $d, (Get-Date -Format o)) | Add-Content -Path $status

  $env:COMPARE_DATASET = $d
  $countCmd = "fprintf('TEMP_COUNT=%d\n', get_dataset_temp_count(getenv('COMPARE_DATASET')));"
  try {
    $out = & $matlab -wait -batch $countCmd 2>&1
  } catch {
    $out = @($_.Exception.Message)
  }

  $m = [regex]::Match(($out -join "`n"), 'TEMP_COUNT=(\d+)')
  if(-not $m.Success){
    ('FAILED count parse: {0} at {1}' -f $d, (Get-Date -Format o)) | Add-Content -Path $status
    continue
  }
  $count = [int]$m.Groups[1].Value
  ('Total temperatures for {0}: {1}' -f $d, $count) | Add-Content -Path $status

  for($start=1; $start -le $count; $start += 10){
    $end = [Math]::Min($start + 9, $count)
    $tag = "b${start}_$end"
    $base = [System.IO.Path]::GetFileNameWithoutExtension($d).ToLower()
    $batchOut = Join-Path $repo ($base + '_paper_vs_residual_peak_cv_comparison_' + $tag + '.txt')

    if(Test-Path $batchOut){
      "  Batch $start-$end already complete, skipping at $(Get-Date -Format o)" | Add-Content -Path $status
      continue
    }

    "  Batch $start-$end start at $(Get-Date -Format o)" | Add-Content -Path $status

    $env:COMPARE_START_IDX = [string]$start
    $env:COMPARE_MAX_TEMPS = '10'
    $env:COMPARE_OUT_TAG = $tag

    $ok = $false
    for($attempt=1; $attempt -le 2; $attempt++){
      "    Attempt $attempt" | Add-Content -Path $status

      # Ensure clean process state before each batch attempt.
      Get-Process -Name MATLAB -ErrorAction SilentlyContinue | Stop-Process -Force -ErrorAction SilentlyContinue

      $logOut = Join-Path $logDir ($base + "_${tag}_attempt${attempt}_out.log")
      $logErr = Join-Path $logDir ($base + "_${tag}_attempt${attempt}_err.log")

      $p = Start-Process -FilePath $matlab -ArgumentList '-wait','-batch','run_paper_vs_residual_peak_cv_compare' -PassThru -NoNewWindow -RedirectStandardOutput $logOut -RedirectStandardError $logErr
      $timedOut = $false
      try {
        Wait-Process -Id $p.Id -Timeout $batchTimeoutSec -ErrorAction Stop
      } catch {
        $timedOut = $true
      }

      if($timedOut){
        Stop-Process -Id $p.Id -Force -ErrorAction SilentlyContinue
        Get-Process -Name MATLAB -ErrorAction SilentlyContinue | Stop-Process -Force -ErrorAction SilentlyContinue
        $code = -999
      } else {
        $p.Refresh()
        $code = $p.ExitCode
        if($null -eq $code -or $code -eq ''){
          if(Test-Path $batchOut){
            $code = 0
          } else {
            $code = -998
          }
        }
      }

      if($code -eq 0){
        $ok = $true
        break
      }

      "    Attempt $attempt failed exit=$code at $(Get-Date -Format o)" | Add-Content -Path $status
      "      Logs: $logOut ; $logErr" | Add-Content -Path $status
    }

    if(-not $ok){
      "  Batch $start-$end FAILED after retries at $(Get-Date -Format o)" | Add-Content -Path $status
      break
    }

    "  Batch $start-$end done at $(Get-Date -Format o)" | Add-Content -Path $status
  }

  Merge-DatasetBatchReports -DatasetFile $d
  "Done dataset: $d at $(Get-Date -Format o)" | Add-Content -Path $status
}

& .\build_all_datasets_compare_summary.ps1
& .\build_all_datasets_lambda_values.ps1
"Completed: $(Get-Date -Format o)" | Add-Content -Path $status
Write-Host 'All datasets complete in 10-temperature batches.'
