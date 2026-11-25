param([switch]$Pozadina)
$ErrorActionPreference = 'Stop'
$projekt = $PSScriptRoot
$python = Join-Path (Split-Path -Parent $projekt) '.venv\Scripts\python.exe'
if (-not (Test-Path -LiteralPath $python)) {
    $python = (Get-Command python -ErrorAction Stop).Source
}
$env:PYTHONUTF8 = '1'
if ($Pozadina) {
    $stdout = Join-Path $projekt ('serija_' + (Get-Date -Format 'yyyyMMdd_HHmmss') + '.log')
    $stderr = $stdout + '.err'
    $proc = Start-Process -FilePath $python -ArgumentList @('-u', 'studija.py', 'run') -WorkingDirectory $projekt -WindowStyle Hidden -RedirectStandardOutput $stdout -RedirectStandardError $stderr -PassThru
    Write-Output "PID $($proc.Id), dnevnik $stdout"
} else {
    Push-Location -LiteralPath $projekt
    try {
        & $python -u studija.py run
        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    } finally {
        Pop-Location
    }
}
