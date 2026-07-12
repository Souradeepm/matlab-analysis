$ErrorActionPreference = 'Stop'
$procs = Get-CimInstance Win32_Process | Where-Object {
  $_.Name -match 'powershell.exe' -and $_.CommandLine -match 'run_all_datasets'
}
foreach($p in $procs){
  Stop-Process -Id $p.ProcessId -Force
}
Write-Output ($procs | Measure-Object).Count
