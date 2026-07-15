$ErrorActionPreference = 'Stop'
Set-Location (Split-Path -Parent $MyInvocation.MyCommand.Path)

$env:BAYES_DRT_MAX_TEMPS = '10'
$env:BAYES_DRT_LAMBDA_MIN_EXP = '-8'
$env:BAYES_DRT_LAMBDA_MAX_EXP = '2'
$env:BAYES_DRT_LAMBDA_COUNT = '5'
$env:BAYES_DRT_DATASET = 'S2022Sap.xlsx'

& 'C:\Progra~1\MATLAB\R2023a\bin\matlab.exe' -wait -batch "disp(pwd); run_bayes_drt_workflow_matlab2011('S2022Sap.xlsx')"
