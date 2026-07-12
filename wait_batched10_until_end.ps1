$ErrorActionPreference = 'Stop'
Set-Location (Split-Path -Parent $MyInvocation.MyCommand.Path)

$procs = Get-CimInstance Win32_Process |
  Where-Object { $_.Name -eq 'powershell.exe' -and $_.CommandLine -match 'run_all_datasets_all_temps_batched10\.ps1' }

if(-not $procs){
  Write-Host 'No active orchestrator process found.'
} else {
  $ids = @($procs.ProcessId)
  Write-Host ('Waiting for orchestrator PIDs: ' + ($ids -join ', '))
  Wait-Process -Id $ids
  Write-Host 'Orchestrator finished.'
}

$statusPath = Join-Path (Get-Location) 'all_datasets_all_temps_batched10_status.txt'
$s = Get-Content $statusPath

"CompletedMarkers=$((($s | Select-String '^Completed:').Count))"
"FailedMarkers=$((($s | Select-String 'FAILED').Count))"
"DoneDatasetLines=$((($s | Select-String '^Done dataset:').Count))"
'---TAIL---'
$s | Select-Object -Last 220
