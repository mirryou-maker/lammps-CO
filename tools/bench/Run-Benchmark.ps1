<#
.SYNOPSIS
Run a LAMMPS input script several times and record the timing breakdown.

.DESCRIPTION
Runs `<Lmp> -in <InputFile> [-var ... ] [ExtraArgs]` `Repeats` times, parses
the "Loop time" / "MPI task timing breakdown" block from each run via
LammpsLog.psm1, and writes:
  - results/logs/<Label>-<timestamp>-run<N>.log   (raw stdout per run)
  - results/<Label>-<timestamp>.json              (full + aggregated results)
  - results/summary.csv                           (one row appended per call)

.EXAMPLE
.\Run-Benchmark.ps1 -Lmp ..\..\build\serial\Release\lmp.exe `
    -InputFile ..\..\lammps-src\examples\melt\in.melt -Label baseline-serial -Repeats 5

.EXAMPLE
.\Run-Benchmark.ps1 -Lmp ..\..\build\omp\Release\lmp.exe `
    -InputFile ..\..\lammps-src\examples\melt\in.melt -Label omp-4t -Repeats 5 `
    -Threads 4 -ExtraArgs '-sf','omp','-pk','omp','4'
#>

[CmdletBinding()]
param(
    [Parameter(Mandatory)] [string]$Lmp,
    [Parameter(Mandatory)] [string]$InputFile,
    [Parameter(Mandatory)] [string]$Label,
    [int]$Repeats = 5,
    [int]$Threads = 0,
    [string[]]$ExtraArgs = @(),
    [string[]]$Vars = @(),
    [string]$ResultsDir,
    [string]$Note = ""
)

$ErrorActionPreference = "Stop"

function Get-MeanStd {
    param([AllowNull()][double[]]$Values)
    $vals = @($Values | Where-Object { $null -ne $_ })
    if ($vals.Count -eq 0) { return [ordered]@{ Mean = $null; Std = $null; N = 0 } }
    $mean = ($vals | Measure-Object -Average).Average
    if ($vals.Count -lt 2) { return [ordered]@{ Mean = $mean; Std = 0.0; N = $vals.Count } }
    $sumSq = ($vals | ForEach-Object { [math]::Pow($_ - $mean, 2) } | Measure-Object -Sum).Sum
    $std = [math]::Sqrt($sumSq / ($vals.Count - 1))
    return [ordered]@{ Mean = $mean; Std = $std; N = $vals.Count }
}

$ScriptDir = $PSScriptRoot
if (-not $ResultsDir) { $ResultsDir = Join-Path $ScriptDir "results" }
$LogsDir = Join-Path $ResultsDir "logs"
New-Item -ItemType Directory -Force -Path $LogsDir | Out-Null

Import-Module (Join-Path $ScriptDir "LammpsLog.psm1") -Force

$Lmp = (Resolve-Path $Lmp).Path
$InputFile = (Resolve-Path $InputFile).Path

$argList = @("-in", $InputFile, "-log", "none")
foreach ($v in $Vars) {
    $parts = $v -split '=', 2
    if ($parts.Count -ne 2) { throw "Invalid -Vars entry '$v', expected NAME=VALUE" }
    $argList += @("-var", $parts[0], $parts[1])
}
$argList += $ExtraArgs

$commit = $null
# Try common LAMMPS clone locations (lammps-src/ or lammps/ beside repo root)
$RepoRoot  = (Resolve-Path (Join-Path $ScriptDir "..\..")).Path
foreach ($candidate in @("lammps-src", "lammps")) {
    $lammpsSrc = Join-Path $RepoRoot $candidate
    if (Test-Path (Join-Path $lammpsSrc ".git")) {
        $commit = (& git -C $lammpsSrc rev-parse --short HEAD 2>$null)
        break
    }
}

$timestamp = Get-Date -Format "yyyyMMdd-HHmmss"

$envBackup = $env:OMP_NUM_THREADS
if ($Threads -gt 0) {
    $env:OMP_NUM_THREADS = "$Threads"
}

$runs = @()
try {
    for ($i = 1; $i -le $Repeats; $i++) {
        Write-Host "[$Label] run $i/$Repeats ..."
        $sw = [System.Diagnostics.Stopwatch]::StartNew()
        $output = & $Lmp @argList 2>&1
        $sw.Stop()

        $logFile = Join-Path $LogsDir "$Label-$timestamp-run${i}.log"
        $output | Out-File -FilePath $logFile -Encoding utf8

        $parsed = $output | ConvertFrom-LammpsLog
        if ($null -eq $parsed.LoopTime) {
            Write-Warning "[$Label] run ${i}: could not find 'Loop time of ...' in output (see $logFile)"
        }

        $runs += [PSCustomObject]@{
            WallTimeSec = $sw.Elapsed.TotalSeconds
            Parsed      = $parsed
        }
    }
}
finally {
    if ($Threads -gt 0) { $env:OMP_NUM_THREADS = $envBackup }
}

$loopStats = Get-MeanStd -Values @($runs | ForEach-Object { $_.Parsed.LoopTime })
$wallStats = Get-MeanStd -Values @($runs | ForEach-Object { $_.WallTimeSec })

$sectionStats = [ordered]@{}
foreach ($name in @('Pair', 'Bond', 'Angle', 'Dihedral', 'Improper', 'Kspace', 'Neigh', 'Comm', 'Output', 'Modify', 'Other')) {
    $avgVals = @($runs | ForEach-Object { $_.Parsed.Sections[$name].AvgTime })
    $pctVals = @($runs | ForEach-Object { $_.Parsed.Sections[$name].PctTotal })
    if (@($avgVals | Where-Object { $null -ne $_ }).Count -gt 0) {
        $sectionStats[$name] = [ordered]@{
            AvgTime  = Get-MeanStd -Values $avgVals
            PctTotal = Get-MeanStd -Values $pctVals
        }
    }
}

$first = $runs[0].Parsed

$resultObj = [ordered]@{
    Label         = $Label
    Timestamp     = $timestamp
    Note          = $Note
    LammpsVersion = $first.Version
    LammpsCommit  = $commit
    LmpPath       = $Lmp
    InputFile     = $InputFile
    ExtraArgs     = $ExtraArgs
    Vars          = $Vars
    ThreadsArg    = $Threads
    Repeats       = $Repeats
    NAtoms        = $first.NAtoms
    NSteps        = $first.NSteps
    NProcs        = $first.NProcs
    OmpThreads    = $first.OmpThreads
    LoopTime      = $loopStats
    WallTime      = $wallStats
    Sections      = $sectionStats
    Runs          = $runs | ForEach-Object {
        [ordered]@{
            WallTimeSec = $_.WallTimeSec
            LoopTime    = $_.Parsed.LoopTime
            Sections    = $_.Parsed.Sections
        }
    }
}

$jsonPath = Join-Path $ResultsDir "$Label-$timestamp.json"
$resultObj | ConvertTo-Json -Depth 10 | Out-File -FilePath $jsonPath -Encoding utf8

$csvRow = [ordered]@{
    Timestamp     = $timestamp
    Label         = $Label
    LammpsVersion = $first.Version
    LammpsCommit  = $commit
    NAtoms        = $first.NAtoms
    NSteps        = $first.NSteps
    NProcs        = $first.NProcs
    ThreadsArg    = $Threads
    OmpThreads    = $first.OmpThreads
    Repeats       = $Repeats
    LoopTimeMean  = $loopStats.Mean
    LoopTimeStd   = $loopStats.Std
    WallTimeMean  = $wallStats.Mean
    WallTimeStd   = $wallStats.Std
}
foreach ($name in @('Pair', 'Bond', 'Angle', 'Dihedral', 'Improper', 'Kspace', 'Neigh', 'Comm', 'Output', 'Modify', 'Other')) {
    $key = "${name}PctMean"
    if ($sectionStats.Contains($name)) {
        $csvRow[$key] = $sectionStats[$name].PctTotal.Mean
    }
    else {
        $csvRow[$key] = $null
    }
}
$csvRow["Note"] = $Note

$csvPath = Join-Path $ResultsDir "summary.csv"
[PSCustomObject]$csvRow | Export-Csv -Path $csvPath -Append -NoTypeInformation -Encoding utf8

Write-Host ""
if ($null -ne $loopStats.Mean) {
    Write-Host ("=== {0}: Loop time = {1:F6} s (std {2:F6}) over {3} run(s), {4} atoms / {5} steps ===" -f `
            $Label, $loopStats.Mean, $loopStats.Std, $loopStats.N, $first.NAtoms, $first.NSteps)
    foreach ($name in $sectionStats.Keys) {
        $s = $sectionStats[$name]
        Write-Host ("  {0,-10} {1,6:F2}%  (avg {2:F6} s)" -f $name, $s.PctTotal.Mean, $s.AvgTime.Mean)
    }
}
else {
    Write-Warning "No successful runs parsed for '$Label' - check logs in $LogsDir"
}
Write-Host ""
Write-Host "Saved: $jsonPath"
Write-Host "Appended: $csvPath"
