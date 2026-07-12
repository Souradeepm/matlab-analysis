$ErrorActionPreference = 'Stop'
Set-Location (Split-Path -Parent $MyInvocation.MyCommand.Path)

$statusPath = Join-Path (Get-Location) 'all_datasets_all_temps_batched10_status.txt'
$s = Get-Content $statusPath

"CompletedMarkers=$((($s | Select-String '^Completed:').Count))"
"FailedMarkers=$((($s | Select-String 'FAILED').Count))"
"DoneDatasetLines=$((($s | Select-String '^Done dataset:').Count))"
'---TAIL---'
$s | Select-Object -Last 120

'---ACTIVE ORCHESTRATOR---'
Get-CimInstance Win32_Process |
  Where-Object { $_.Name -eq 'powershell.exe' -and $_.CommandLine -match 'run_all_datasets_all_temps_batched10\.ps1' } |
  Select-Object ProcessId, Name, CommandLine |
  Format-List

'---ACTIVE MATLAB---'
Get-Process -Name matlab, MATLAB -ErrorAction SilentlyContinue |
  Select-Object Name, Id, CPU, StartTime, Responding |
  Format-Table -AutoSize
