param(
    [int]   $OmpThreads = 4,
    [int]   $Steps      = 100000,
    [switch]$DryRun,
    [string]$LogDir     = "tools\fep\logs_opls_long"
)

Set-Location (Split-Path (Split-Path $PSScriptRoot -Parent) -Parent)
$lmp = "build\omp-fep\Release\lmp.exe"
New-Item -ItemType Directory -Force -Path $LogDir | Out-Null

$lambdas = @(0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,
             0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.00)

$solvents = @(
    @{name='ec';  infile='tools\fep\in.fep_ec_pc'; data='tools\fep\data\ec_li_aa.data'},
    @{name='pc';  infile='tools\fep\in.fep_ec_pc'; data='tools\fep\data\pc_li_aa.data'},
    @{name='dme'; infile='tools\fep\in.fep_dme';   data='tools\fep\data\dme_li_aa.data'}
)

$total   = $lambdas.Count * $solvents.Count
$done    = 0
$t_start = Get-Date

foreach ($sol in $solvents) {
    $s = $sol['name']
    foreach ($lam in $lambdas) {
        $lamStr = "{0:F2}" -f $lam
        $tag    = "${s}_lam${lamStr}"
        $stdout = "$LogDir\fep_${tag}.stdout"
        $stderr = "$LogDir\fep_${tag}.stderr"

        # Skip if already completed (has Total wall time)
        if (Test-Path $stdout) {
            $content = Get-Content $stdout -Raw -ErrorAction SilentlyContinue
            if ($content -match "Total wall time") { $done++; continue }
        }

        $args_ = @(
            "-sf","omp","-pk","omp",$OmpThreads,
            "-var","DATAFILE",$sol['data'],
            "-var","LAMBDA",$lamStr,
            "-var","NSTEPS",$Steps,
            "-in",$sol['infile'],
            "-log","none"
        )

        if ($DryRun) {
            Write-Host "DRY $tag  -> $stdout"
        } else {
            $t0 = Get-Date
            & $lmp @args_ 2>"$stderr" | Out-File -FilePath "$stdout" -Encoding UTF8
            $el  = [math]::Round(((Get-Date)-$t0).TotalSeconds,1)
            $done++
            $pct = [math]::Round(100*$done/$total)
            Write-Host ("[{0,3}%] {1,-20}  {2,5}s" -f $pct,$tag,$el)
        }
    }
}

$tot = [math]::Round(((Get-Date)-$t_start).TotalSeconds/60,1)
Write-Host "Completed $done/$total windows in ${tot} min."
