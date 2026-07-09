$ErrorActionPreference = 'Stop'
$root = 'c:\Users\mitra\Downloads\matlab analysis'
$file = Join-Path $root 'lambda_perturbation_result.xlsx'

$excel = New-Object -ComObject Excel.Application
$excel.Visible = $false
$excel.DisplayAlerts = $false

if (-not (Test-Path $file)) {
    throw "File not found: $file"
}

$wb = $excel.Workbooks.Open($file, $null, $true)
try {
    $rows = @()
    foreach ($ws in $wb.Worksheets) {
        $sheetName = [string]$ws.Name
        if ($sheetName -notlike 'LAMBDA_*') { continue }

        $vals = $ws.UsedRange.Value2
        if ($null -eq $vals) { continue }

        $nRows = $vals.GetLength(0)
        $nCols = $vals.GetLength(1)
        if ($nRows -lt 2 -or $nCols -lt 5) { continue }

        $delta = @()
        $lambdaRef = $null
        for ($r = 2; $r -le $nRows; $r++) {
            $factor = $vals[$r,1]
            $lambda = $vals[$r,2]
            $d = $vals[$r,5]
            if ($null -ne $d) { $delta += [double]$d }
            if ($null -ne $factor -and [double]$factor -eq 1.0) {
                $lambdaRef = [double]$lambda
            }
        }

        if ($delta.Count -eq 0) { continue }

        $tempText = $sheetName -replace '^LAMBDA_', ''
        $tempText = $tempText -replace 'p', '.'
        $tempText = $tempText -replace 'K$', ''
        $tempVal = $null
        [double]::TryParse($tempText, [ref]$tempVal) | Out-Null

        $rows += [pscustomobject]@{
            Temperature_K = $tempVal
            Sheet = $sheetName
            Lambda_ref = $lambdaRef
            Min_delta_pct = ($delta | Measure-Object -Minimum).Minimum
            Max_delta_pct = ($delta | Measure-Object -Maximum).Maximum
            MaxAbs_delta_pct = ($delta | ForEach-Object { [math]::Abs($_) } | Measure-Object -Maximum).Maximum
        }
    }

    $rows | Sort-Object Temperature_K | ConvertTo-Csv -NoTypeInformation
}
finally {
    $wb.Close($false)
    $excel.Quit()
    [System.Runtime.InteropServices.Marshal]::ReleaseComObject($wb) | Out-Null
    [System.Runtime.InteropServices.Marshal]::ReleaseComObject($excel) | Out-Null
}
