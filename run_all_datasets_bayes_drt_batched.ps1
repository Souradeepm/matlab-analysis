$ErrorActionPreference = 'Stop'

$repo = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $repo

# ---- Configuration -------------------------------------------------------
$env:BAYES_DRT_MAX_TEMPS     = '10'
$env:BAYES_DRT_LAMBDA_MIN_EXP = '-8'
$env:BAYES_DRT_LAMBDA_MAX_EXP = '2'
$env:BAYES_DRT_LAMBDA_COUNT  = '21'

$batchSize      = 10
$batchTimeoutSec = 480   # seconds per batch attempt; override via BAYES_DRT_BATCH_TIMEOUT_SEC
$tmp = 0.0
if([double]::TryParse($env:BAYES_DRT_BATCH_TIMEOUT_SEC,[ref]$tmp) -and $tmp -ge 60){
  $batchTimeoutSec = [int][Math]::Round($tmp)
}

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

$matlab  = 'C:\Progra~1\MATLAB\R2023a\bin\matlab.exe'
$status  = Join-Path $repo 'all_datasets_bayes_drt_batched_status.txt'
$logDir  = Join-Path $repo 'bayes_drt_batch_logs'
$summaryFile = Join-Path $repo 'all_datasets_bayes_drt_all_temps_summary.txt'

"Start: $(Get-Date -Format o)" | Set-Content -Path $status -Encoding ascii
if(-not (Test-Path $logDir)){ New-Item -ItemType Directory -Path $logDir | Out-Null }

# ---- Merge function -------------------------------------------------------
# Reads all batch fragment files for a dataset, merges per-temp CSV rows,
# writes the merged per-dataset report, and returns a summary PSCustomObject.
function Merge-BayesDrtBatchReports {
  param([string]$DatasetFile)

  $base = [System.IO.Path]::GetFileNameWithoutExtension($DatasetFile).ToLower()
  $pattern = $base + '_bayes_drt_matlab2011_b*_*.txt'
  $files = Get-ChildItem -Path $repo -Filter $pattern -File -ErrorAction SilentlyContinue |
    Sort-Object {
      $m = [regex]::Match($_.BaseName, '_b(\d+)_')
      if($m.Success){ [int]$m.Groups[1].Value } else { 999999 }
    }

  if($files.Count -eq 0){ return $null }

  $rows = @()
  foreach($f in $files){
    $lines = Get-Content $f.FullName
    $headerIdx = -1
    for($i = 0; $i -lt $lines.Count; $i++){
      if($lines[$i] -match '^TemperatureK,SelectedLambda'){
        $headerIdx = $i; break
      }
    }
    if($headerIdx -lt 0){ continue }
    for($j = $headerIdx + 1; $j -lt $lines.Count; $j++){
      $line = $lines[$j].Trim()
      if([string]::IsNullOrWhiteSpace($line)){ continue }
      $p = $line.Split(',')
      if($p.Count -lt 7){ continue }
      $t = 0.0; $lam = 0.0; $rc = 0.0; $ic = 0.0; $tc = 0.0; $mr = 0.0; $pk = 0
      if(-not [double]::TryParse($p[0],[ref]$t)){ continue }
      if(-not [double]::TryParse($p[1],[ref]$lam)){ continue }
      if(-not [double]::TryParse($p[2],[ref]$rc)){ continue }
      if(-not [double]::TryParse($p[3],[ref]$ic)){ continue }
      if(-not [double]::TryParse($p[4],[ref]$tc)){ continue }
      if(-not [double]::TryParse($p[5],[ref]$mr)){ continue }
      if(-not [int]::TryParse($p[6],[ref]$pk)){ continue }
      $rows += [PSCustomObject]@{
        TemperatureK    = $t
        SelectedLambda  = $lam
        RealCV          = $rc
        ImagCV          = $ic
        TotalCV         = $tc
        MeanAbsResidual = $mr
        PeakCount       = $pk
      }
    }
  }

  if($rows.Count -eq 0){ return $null }
  $rows = $rows | Sort-Object TemperatureK

  $meanLambda  = ($rows | Measure-Object -Property SelectedLambda  -Average).Average
  $meanRealCV  = ($rows | Measure-Object -Property RealCV          -Average).Average
  $meanImagCV  = ($rows | Measure-Object -Property ImagCV          -Average).Average
  $meanTotalCV = ($rows | Measure-Object -Property TotalCV         -Average).Average
  $meanResid   = ($rows | Measure-Object -Property MeanAbsResidual -Average).Average
  $meanPeaks   = ($rows | Measure-Object -Property PeakCount       -Average).Average
  $maxPeaks    = ($rows | Measure-Object -Property PeakCount       -Maximum).Maximum

  $out = Join-Path $repo ($base + '_bayes_drt_matlab2011_all_temps.txt')
  $o   = @()
  $o  += 'Bayes-DRT MATLAB 2011 workflow (Re/Im CV) — merged all temperatures'
  $o  += "Dataset: $DatasetFile"
  $o  += "Generated: $(Get-Date -Format 'dd-MMM-yyyy HH:mm:ss')"
  $o  += "Total temperatures: $($rows.Count)"
  $o  += ''
  $o  += 'Aggregate summary'
  $o  += ('Mean selected lambda: {0:E8}'  -f $meanLambda)
  $o  += ('Mean real-part CV: {0:E8}'     -f $meanRealCV)
  $o  += ('Mean imag-part CV: {0:E8}'     -f $meanImagCV)
  $o  += ('Mean total CV: {0:E8}'         -f $meanTotalCV)
  $o  += ('Mean absolute residual: {0:E8}' -f $meanResid)
  $o  += ('Mean peak count: {0:F4}'       -f $meanPeaks)
  $o  += ('Max peak count: {0}'           -f $maxPeaks)
  $o  += ''
  $o  += 'Per-temperature detail'
  $o  += 'TemperatureK,SelectedLambda,RealCV,ImagCV,TotalCV,MeanAbsResidual,PeakCount'
  foreach($r in $rows){
    $o += ('{0:F6},{1:E8},{2:E8},{3:E8},{4:E8},{5:E8},{6}' -f `
      $r.TemperatureK,$r.SelectedLambda,$r.RealCV,$r.ImagCV,$r.TotalCV,$r.MeanAbsResidual,$r.PeakCount)
  }
  Set-Content -Path $out -Value $o -Encoding ascii

  return [PSCustomObject]@{
    Dataset         = $base
    Temperatures    = $rows.Count
    MeanLambda      = $meanLambda
    MeanRealCV      = $meanRealCV
    MeanImagCV      = $meanImagCV
    MeanTotalCV     = $meanTotalCV
    MeanResidual    = $meanResid
    MeanPeakCount   = $meanPeaks
    MaxPeakCount    = $maxPeaks
    ReportFile      = $out
  }
}

# ---- Main loop -----------------------------------------------------------
$summaryRows = @()

foreach($d in $datasets){

  "Running dataset: $d at $(Get-Date -Format o)" | Add-Content -Path $status

  $env:BAYES_DRT_DATASET = $d

  # Get temperature count via MATLAB
  $countCmd = "fprintf('TEMP_COUNT=%d\n', get_dataset_temp_count(getenv('BAYES_DRT_DATASET')));"
  try { $out = & $matlab -wait -batch $countCmd 2>&1 } catch { $out = @($_.Exception.Message) }
  $m = [regex]::Match(($out -join "`n"), 'TEMP_COUNT=(\d+)')
  if(-not $m.Success){
    "FAILED count parse for $d" | Add-Content -Path $status
    continue
  }
  $count = [int]$m.Groups[1].Value
  "Total temperatures for ${d}: $count" | Add-Content -Path $status

  for($start = 1; $start -le $count; $start += $batchSize){
    $end = [Math]::Min($start + $batchSize - 1, $count)
    $tag = "b${start}_$end"
    $base = [System.IO.Path]::GetFileNameWithoutExtension($d).ToLower()
    $batchOut = Join-Path $repo ($base + '_bayes_drt_matlab2011_' + $tag + '.txt')

    # Resume-safe: skip if batch output already exists
    if(Test-Path $batchOut){
      "  Batch $start-$end already complete, skipping" | Add-Content -Path $status
      continue
    }

    "  Batch $start-$end start at $(Get-Date -Format o)" | Add-Content -Path $status

    $env:BAYES_DRT_START_IDX = [string]$start
    $env:BAYES_DRT_MAX_TEMPS  = [string]$batchSize
    $env:BAYES_DRT_OUT_TAG    = $tag

    $ok = $false
    for($attempt = 1; $attempt -le 2; $attempt++){
      "    Attempt $attempt" | Add-Content -Path $status

      Get-Process -Name MATLAB -ErrorAction SilentlyContinue | Stop-Process -Force -ErrorAction SilentlyContinue

      $outLog = Join-Path $logDir ($base + "_${tag}_attempt${attempt}_out.log")
      $errLog = Join-Path $logDir ($base + "_${tag}_attempt${attempt}_err.log")

      $p = Start-Process -FilePath $matlab `
            -ArgumentList '-wait','-batch','run_bayes_drt_workflow_matlab2011' `
            -PassThru -NoNewWindow `
            -RedirectStandardOutput $outLog -RedirectStandardError $errLog

      $timedOut = $false
      try { Wait-Process -Id $p.Id -Timeout $batchTimeoutSec -ErrorAction Stop }
      catch { $timedOut = $true }

      if($timedOut){
        Stop-Process -Id $p.Id -Force -ErrorAction SilentlyContinue
        Get-Process -Name MATLAB -ErrorAction SilentlyContinue | Stop-Process -Force -ErrorAction SilentlyContinue
        $code = -999
      } else {
        $p.Refresh(); $code = $p.ExitCode
        if($null -eq $code -or $code -eq ''){
          $code = if(Test-Path $batchOut){ 0 } else { -998 }
        }
      }

      if($code -eq 0 -and (Test-Path $batchOut)){
        $ok = $true
        "    Attempt $attempt success" | Add-Content -Path $status
        break
      }
      "    Attempt $attempt failed exit=$code" | Add-Content -Path $status
    }

    if(-not $ok){
      "  Batch $start-$end FAILED after retries at $(Get-Date -Format o)" | Add-Content -Path $status
      break
    }
    "  Batch $start-$end done at $(Get-Date -Format o)" | Add-Content -Path $status
  }

  # Merge all batch fragments into a single per-dataset report
  $merged = Merge-BayesDrtBatchReports -DatasetFile $d
  if($null -ne $merged){
    $summaryRows += $merged
    "Done dataset: $d at $(Get-Date -Format o)" | Add-Content -Path $status
  } else {
    "WARNING: merge produced no rows for $d" | Add-Content -Path $status
  }
}

# ---- Write aggregate summary --------------------------------------------
$summaryLines = @()
$summaryLines += 'Dataset,Temperatures,MeanSelectedLambda,MeanRealCV,MeanImagCV,MeanTotalCV,MeanResidual,MeanPeakCount,MaxPeakCount,ReportFile'
foreach($row in $summaryRows){
  $summaryLines += ('{0},{1},{2:E8},{3:E8},{4:E8},{5:E8},{6:E8},{7:F4},{8},{9}' -f `
    $row.Dataset,$row.Temperatures,$row.MeanLambda,$row.MeanRealCV,$row.MeanImagCV,`
    $row.MeanTotalCV,$row.MeanResidual,$row.MeanPeakCount,$row.MaxPeakCount,$row.ReportFile)
}
$summaryLines | Set-Content -Path $summaryFile -Encoding ascii

"Completed: $(Get-Date -Format o)" | Add-Content -Path $status
Write-Host 'Bayes-DRT all-temperatures batched run complete.'
