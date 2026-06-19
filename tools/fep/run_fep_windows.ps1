# ─────────────────────────────────────────────────────────────────────────────
# run_fep_windows.ps1  — Run 20-window FEP for Li+ desolvation
# Usage:
#   .\tools\fep\run_fep_windows.ps1 [-OmpThreads 4] [-Steps 1000000] [-DryRun]
#
# From: D:\Claude-Code-R\LAMMPS-CO\  (project root)
# ─────────────────────────────────────────────────────────────────────────────

param(
    [int]$OmpThreads = 4,
    [int]$Steps      = 500000,    # production steps per window (1 ns = 500k x 2fs)
    [switch]$DryRun,
    [string]$LogDir  = "tools\fep\logs"
)

$lmp     = "build\omp-fep\Release\lmp.exe"
$infile  = "tools\fep\in.fep_lj"
$lambdas = @(0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
             0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00)

New-Item -ItemType Directory -Force $LogDir | Out-Null

$total_start = Get-Date
Write-Host "FEP 20-window run: $($lambdas.Count) windows, $OmpThreads OMP threads, $Steps steps/window"
Write-Host "Output: $LogDir"
Write-Host ("-"*65)

$results = @()

foreach ($lam in $lambdas) {
    $label   = $lam.ToString("0.00")
    $logfile = "$LogDir\fep_lambda${label}.log"
    $t0      = Get-Date

    Write-Host -NoNewline "  lambda=$label  "

    if ($DryRun) {
        Write-Host "(dry run - skipped)"
        continue
    }

    # Build LAMMPS argument list
    $args_list = @(
        "-sf", "omp", "-pk", "omp", "$OmpThreads",
        "-var", "LAMBDA", "$lam",
        "-var", "NSTEPS", "$Steps",
        "-in",  $infile,
        "-log", $logfile
    )

    $proc = Start-Process -FilePath $lmp -ArgumentList $args_list `
                          -NoNewWindow -Wait -PassThru `
                          -RedirectStandardOutput "$LogDir\fep_lambda${label}.stdout" `
                          -RedirectStandardError  "$LogDir\fep_lambda${label}.stderr"

    $elapsed = [math]::Round(((Get-Date) - $t0).TotalSeconds, 1)

    if ($proc.ExitCode -eq 0) {
        # Quick parse: count thermo lines (production data)
        $n_frames = (Get-Content "$LogDir\fep_lambda${label}.stdout" |
                     Where-Object { $_ -match "^\s+[0-9]" }).Count
        Write-Host "OK  ($elapsed s, $n_frames frames)"
        $results += [PSCustomObject]@{Lambda=$lam; Status="OK"; Time=$elapsed; Frames=$n_frames}
    } else {
        $err = (Get-Content "$LogDir\fep_lambda${label}.stderr" | Select-Object -Last 3) -join " | "
        Write-Host "FAIL  ($err)"
        $results += [PSCustomObject]@{Lambda=$lam; Status="FAIL"; Time=$elapsed; Frames=0}
    }
}

$total_elapsed = [math]::Round(((Get-Date) - $total_start).TotalSeconds / 60, 1)
Write-Host ("-"*65)
Write-Host "Total: $total_elapsed min  |  $($results | Where-Object Status -eq 'OK' | Measure-Object).Count / $($lambdas.Count) windows OK"

# Save summary
$results | Export-Csv "$LogDir\run_summary.csv" -NoTypeInformation
Write-Host "Summary → $LogDir\run_summary.csv"
