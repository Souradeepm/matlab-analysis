$ErrorActionPreference = 'Stop'
$root = 'c:\Users\mitra\Downloads\matlab analysis'

$excel = New-Object -ComObject Excel.Application
$excel.Visible = $false
$excel.DisplayAlerts = $false

function Read-Summary {
    param(
        [string]$File,
        [string]$Sheet
    )

    if (-not (Test-Path $File)) {
        return @()
    }

    $wb = $excel.Workbooks.Open($File, $null, $true)
    try {
        $ws = $wb.Worksheets.Item($Sheet)
        $vals = $ws.UsedRange.Value2
        if ($null -eq $vals) {
            return @()
        }

        $rows = $vals.GetLength(0)
        $cols = $vals.GetLength(1)
        if ($rows -lt 2) {
            return @()
        }

        $headers = @()
        for ($c = 1; $c -le $cols; $c++) {
            $headers += [string]$vals[1, $c]
        }

        $items = @()
        for ($r = 2; $r -le $rows; $r++) {
            $obj = [ordered]@{}
            for ($c = 1; $c -le $cols; $c++) {
                $h = $headers[$c - 1]
                if ([string]::IsNullOrWhiteSpace($h)) {
                    $h = "col_$c"
                }
                $obj[$h] = $vals[$r, $c]
            }
            $items += New-Object psobject -Property $obj
        }

        return $items
    }
    finally {
        $wb.Close($false)
    }
}

$kk = Read-Summary -File (Join-Path $root 'kk_result.xlsx') -Sheet "summary"
$sens = Read-Summary -File (Join-Path $root 'sensitivity_result.xlsx') -Sheet "summary"
$lam = Read-Summary -File (Join-Path $root 'lambda_perturbation_result.xlsx') -Sheet "summary"

$temps = @($kk.Temperature_K + $sens.Temperature_K + $lam.Temperature_K | Where-Object { $_ -ne $null } | Sort-Object -Unique)
$rows = @()

foreach ($t in $temps) {
    $k = $kk | Where-Object { $_.Temperature_K -eq $t } | Select-Object -First 1
    $s = $sens | Where-Object { $_.Temperature_K -eq $t } | Select-Object -First 1
    $l = $lam | Where-Object { $_.Temperature_K -eq $t } | Select-Object -First 1

    $rows += [pscustomobject]@{
        Temperature_K = $t
        Lambda = if ($l) { $l.Lambda_ref } elseif ($k) { $k.Lambda } else { $null }
        KK_Total_pct = if ($k) { $k.KK_total_pct } else { $null }
        Mean_CV_pct = if ($s) { $s.Mean_CV_pct } else { $null }
        Max_CV_pct = if ($s) { $s.Max_CV_pct } else { $null }
        Lambda_MinDelta_pct = if ($l) { $l.Min_delta_pct } else { $null }
        Lambda_MaxDelta_pct = if ($l) { $l.Max_delta_pct } else { $null }
        Lambda_MaxAbsDelta_pct = if ($l) { $l.Max_abs_delta_pct } else { $null }
    }
}

$rows | Sort-Object Temperature_K | ConvertTo-Csv -NoTypeInformation

$excel.Quit()
[System.Runtime.InteropServices.Marshal]::ReleaseComObject($excel) | Out-Null
