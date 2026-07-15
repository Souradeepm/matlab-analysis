# run_all_datasets_three_methods_batched.ps1
# Orchestrates batch processing of all three lambda selection methods
# across all 11 datasets in batches of 10 temperatures per MATLAB process

$ErrorActionPreference = 'Stop'

$repo = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $repo

# ---- Configuration -------------------------------------------------------
$batchSize        = 10
$batchTimeoutSec  = 480   # 8 minutes per batch
$maxRetries       = 2

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
$status  = Join-Path $repo 'all_datasets_three_methods_batched_status.txt'
$logDir  = Join-Path $repo 'three_methods_batch_logs'

"Start: $(Get-Date -Format o)" | Set-Content -Path $status -Encoding ascii
if(-not (Test-Path $logDir)){ New-Item -ItemType Directory -Path $logDir | Out-Null }

# ---- Helper: Count temperatures in dataset ----
function Get-TemperatureCount {
  param([string]$DatasetFile)
  try {
    $excel = New-Object -ComObject Excel.Application
    $excel.Visible = $false
    $workbook = $excel.Workbooks.Open((Resolve-Path $DatasetFile).Path)
    $sheet = $workbook.Sheets(1)
    $usedRange = $sheet.UsedRange
    $n_cols = $usedRange.Columns.Count
    $workbook.Close()
    $excel.Quit()
    return [Math]::Floor($n_cols / 3)
  } catch {
    Write-Host "Error counting temps in $DatasetFile : $_"
    return 0
  }
}

# ---- Helper: Process batch for all three methods ----
function Process-BatchAllMethods {
  param([string]$DatasetFile, [int]$StartIdx, [int]$EndIdx, [int]$Attempt)
  
  $base = [System.IO.Path]::GetFileNameWithoutExtension($DatasetFile).ToLower()
  $batch = [Math]::Ceiling($StartIdx / 10)
  $logFile = Join-Path $logDir "$($base)_b$($batch)_attempt$($Attempt).log"
  
  Write-Host "  Batch $batch (temps $StartIdx-$EndIdx, attempt $Attempt)..."
  
  $matlab_cmd = @"
addpath('$repo');
env_vars = struct('MATLAB_ANALYSIS_REPO_ROOT', '$repo');
for fname = fieldnames(env_vars)'
  setenv(fname{1}, env_vars.(fname{1}));
end
try
  run_all_three_methods_batch('$DatasetFile', $StartIdx, $EndIdx);
catch ME
  fprintf(2, 'ERROR: %s\n', ME.message);
  exit(1);
end
exit(0);
"@
  
  $startTime = Get-Date
  $proc = Start-Process -FilePath $matlab -ArgumentList "-batch", $matlab_cmd `
    -RedirectStandardOutput $logFile -PassThru -NoNewWindow
  
  $waited = $proc | Wait-Process -Timeout $batchTimeoutSec -PassThru -ErrorAction SilentlyContinue
  
  if (-not $waited) {
    Write-Host "    ✗ Timeout after ${batchTimeoutSec}s. Killing process..."
    Stop-Process -InputObject $proc -Force -ErrorAction SilentlyContinue
    return $false
  }
  
  $elapsed = ((Get-Date) - $startTime).TotalSeconds
  
  if ($proc.ExitCode -eq 0) {
    Write-Host "    ✓ Success (${elapsed}s)"
    return $true
  } else {
    Write-Host "    ✗ Failed with exit code $($proc.ExitCode)"
    Get-Content $logFile -ErrorAction SilentlyContinue | Select-Object -Last 10 | ForEach-Object {
      Write-Host "      $_"
    }
    return $false
  }
}

# ---- Main batch processing loop ----
$totalDatasets = $datasets.Count
$totalBatches = 0
$successBatches = 0
$failedBatches = 0

foreach ($dataset in $datasets) {
  $base = [System.IO.Path]::GetFileNameWithoutExtension($dataset).ToLower()
  
  if (-not (Test-Path $dataset)) {
    Write-Host "✗ Dataset not found: $dataset"
    Add-Content -Path $status -Encoding ascii -Value "SKIP: $base (file not found)"
    continue
  }
  
  Write-Host ""
  Write-Host "Processing $base..."
  $n_temps = Get-TemperatureCount $dataset
  Write-Host "  Temperatures: $n_temps"
  
  if ($n_temps -eq 0) {
    Write-Host "  ✗ Unable to determine temperature count"
    Add-Content -Path $status -Encoding ascii -Value "FAIL: $base (unable to count temperatures)"
    continue
  }
  
  # Process batches
  $start = 1
  while ($start -le $n_temps) {
    $end = [Math]::Min($start + $batchSize - 1, $n_temps)
    $totalBatches += 1
    
    $success = $false
    for ($attempt = 1; $attempt -le $maxRetries; $attempt++) {
      if (Process-BatchAllMethods $dataset $start $end $attempt) {
        $success = $true
        $successBatches += 1
        break
      }
    }
    
    if (-not $success) {
      Write-Host "  ✗ Batch failed after $maxRetries attempts"
      $failedBatches += 1
      Add-Content -Path $status -Encoding ascii -Value "FAIL: $base batch $(($start-1)/10 + 1)"
    }
    
    $start = $end + 1
  }
  
  $n_batches = [Math]::Ceiling($n_temps / $batchSize)
  Add-Content -Path $status -Encoding ascii -Value "OK: $base ($n_temps temps in $n_batches batches)"
}

# ---- Summary ----
Write-Host ""
Write-Host "============================================"
Write-Host "Batch Processing Complete"
Write-Host "============================================"
Write-Host "Total batches: $totalBatches"
Write-Host "Successful:    $successBatches"
Write-Host "Failed:        $failedBatches"
Write-Host "End: $(Get-Date -Format o)"
Write-Host ""

Add-Content -Path $status -Encoding ascii -Value ""
Add-Content -Path $status -Encoding ascii -Value "Total batches: $totalBatches"
Add-Content -Path $status -Encoding ascii -Value "Successful: $successBatches"
Add-Content -Path $status -Encoding ascii -Value "Failed: $failedBatches"
Add-Content -Path $status -Encoding ascii -Value "End: $(Get-Date -Format o)"

if ($failedBatches -gt 0) {
  Write-Host "⚠ Some batches failed. Check log files in: $logDir"
  exit(1)
} else {
  Write-Host "✓ All batches completed successfully"
  exit(0)
}
