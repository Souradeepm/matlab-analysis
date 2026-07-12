$ErrorActionPreference = 'SilentlyContinue'
$procs = Get-CimInstance Win32_Process | Where-Object {
  $_.Name -match 'powershell.exe' -and (
    $_.CommandLine -match 'run_all_datasets_all_temps' -or
    $_.CommandLine -match 'run_all_datasets_full_wait'
  )
}
foreach($p in $procs){ Stop-Process -Id $p.ProcessId -Force }
Write-Output ("KilledLoopProcs=" + ($procs | Measure-Object).Count)
