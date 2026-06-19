<#
.SYNOPSIS
Check a LAMMPS input script for accelerator-package (OMP/OPT/INTEL/GPU/KOKKOS) variants
that exist for the styles it uses but are not currently enabled.

.DESCRIPTION
Scans lammps-src/src/**/*.h for "<Category>Style(name,Class)" registrations to build an
index of which base pair/bond/angle/dihedral/improper/kspace/fix/compute styles have
accelerated (/omp, /opt, /intel, /gpu, /kk[/host|/device]) variants. Then parses the given
LAMMPS input script for pair_style/bond_style/angle_style/dihedral_style/improper_style/
kspace_style/fix/compute commands (including hybrid/hybrid-overlay sub-styles) and any
'suffix' command already present, and reports which used styles have accelerated variants
that are not currently active, with a suggested '-sf <name>' / 'suffix <name>' to enable
them.

This is a static-analysis / informational tool only - it does not modify the input script
or run LAMMPS.

.PARAMETER InputFile
Path to the LAMMPS input script to check (e.g. lammps-src/examples/melt/in.melt).

.PARAMETER BuildDir
Name of the build directory under LAMMPS-CO\build\ whose CMakeCache.txt is checked to
report whether PKG_OPENMP/PKG_OPT/PKG_INTEL/PKG_GPU/PKG_KOKKOS are currently enabled.
Default: "serial".

.PARAMETER Json
Optional path to also write the full report as JSON.

.EXAMPLE
.\Check-AccelCoverage.ps1 -InputFile ..\..\lammps-src\examples\melt\in.melt

.EXAMPLE
.\Check-AccelCoverage.ps1 -InputFile in.foo -BuildDir omp -Json results\foo.json
#>

[CmdletBinding()]
param(
    [Parameter(Mandatory)] [string]$InputFile,
    [string]$BuildDir = "serial",
    [string]$Json
)

$ErrorActionPreference = "Stop"

$ScriptDir = $PSScriptRoot
$RepoRoot  = (Resolve-Path (Join-Path $ScriptDir "..\..")).Path
$SrcDir    = Join-Path $RepoRoot "lammps-src\src"

$Categories = @('Pair', 'Bond', 'Angle', 'Dihedral', 'Improper', 'KSpace', 'Fix', 'Compute')
$AccelDirs  = @{ OPENMP = 'omp'; OPT = 'opt'; INTEL = 'intel'; GPU = 'gpu'; KOKKOS = 'kk' }
$VariantKeys = @('omp', 'opt', 'intel', 'gpu', 'kk')
$PkgNames = [ordered]@{ omp = 'OPENMP'; opt = 'OPT'; intel = 'INTEL'; gpu = 'GPU'; kk = 'KOKKOS' }
$DisplayLabel = @{
    Pair = 'pair_style'; Bond = 'bond_style'; Angle = 'angle_style'; Dihedral = 'dihedral_style'
    Improper = 'improper_style'; KSpace = 'kspace_style'; Fix = 'fix'; Compute = 'compute'
}
$ExtraArgsHint = @{
    omp   = "-pk omp <Nthreads>  (e.g. -sf omp -pk omp 4)"
    opt   = "(no extra -pk needed)"
    intel = "-pk intel 0  (CPU mode)"
    gpu   = "-pk gpu 1  (1 GPU per node)"
    kk    = "-k on t <Nthreads>  (e.g. -k on t 4 -sf kk)"
}

# ---------------------------------------------------------------------------
# 1. Build the style index by scanning lammps-src/src/**/*.h
# ---------------------------------------------------------------------------
function Build-StyleIndex {
    param(
        [string]$SrcDir,
        [string[]]$Categories,
        [hashtable]$AccelDirs
    )

    $pattern = '^\s*(' + ($Categories -join '|') + ')Style\(\s*([^,]+?)\s*,'
    $srcFull = (Resolve-Path $SrcDir).Path.TrimEnd('\', '/')
    $files = Get-ChildItem -Path $srcFull -Recurse -Filter *.h -File
    $hits = $files | Select-String -Pattern $pattern

    $index = @{}
    foreach ($cat in $Categories) { $index[$cat] = @{} }

    foreach ($h in $hits) {
        $m = $h.Matches[0]
        $cat = $m.Groups[1].Value
        $fullName = $m.Groups[2].Value.Trim()
        if ($fullName -eq 'DEPRECATED') { continue }

        $rel = $h.Path.Substring($srcFull.Length).TrimStart('\', '/')
        $topDir = ($rel -split '[\\/]')[0]
        $pkg = 'standard'
        if ($AccelDirs.ContainsKey($topDir)) { $pkg = $AccelDirs[$topDir] }

        $baseName = $fullName
        if ($fullName -match '^(.*)/(kk/host|kk/device|kk|omp|opt|intel|gpu)$') {
            $baseName = $Matches[1]
        }

        if (-not $index[$cat].ContainsKey($baseName)) {
            $index[$cat][$baseName] = @{ standard = $false; omp = $false; opt = $false; intel = $false; gpu = $false; kk = $false }
        }
        $index[$cat][$baseName][$pkg] = $true
    }

    Write-Host "Indexed $($hits.Count) style registrations from $($files.Count) headers under $SrcDir"
    return $index
}

# ---------------------------------------------------------------------------
# 2. Helpers for parsing the input script
# ---------------------------------------------------------------------------
function Get-LogicalLines {
    param([string[]]$Lines)
    $result = @()
    $buffer = ""
    foreach ($raw in $Lines) {
        $line = $raw
        $hashIdx = $line.IndexOf('#')
        if ($hashIdx -ge 0) { $line = $line.Substring(0, $hashIdx) }
        $line = $line.TrimEnd()
        if ($line.EndsWith('&')) {
            $buffer += $line.Substring(0, $line.Length - 1) + ' '
            continue
        }
        $buffer += $line
        $trimmed = $buffer.Trim()
        if ($trimmed.Length -gt 0) { $result += $trimmed }
        $buffer = ""
    }
    return $result
}

function Get-BaseAndSuffix {
    param([string]$Name)
    if ($Name -match '^(.*)/(kk/host|kk/device|kk|omp|opt|intel|gpu)$') {
        $suffix = $Matches[2]
        if ($suffix.StartsWith('kk')) { $suffix = 'kk' }
        return @{ Base = $Matches[1]; Suffix = $suffix }
    }
    return @{ Base = $Name; Suffix = $null }
}

# ---------------------------------------------------------------------------
# 3. Build index, parse input script
# ---------------------------------------------------------------------------
$Index = Build-StyleIndex -SrcDir $SrcDir -Categories $Categories -AccelDirs $AccelDirs

$InputPath = (Resolve-Path $InputFile).Path
$logical = Get-LogicalLines -Lines (Get-Content -Path $InputPath)

$hybridStyles = @('hybrid', 'hybrid/overlay', 'hybrid/scaled')
$catMap = @{
    'pair_style' = 'Pair'; 'bond_style' = 'Bond'; 'angle_style' = 'Angle'
    'dihedral_style' = 'Dihedral'; 'improper_style' = 'Improper'; 'kspace_style' = 'KSpace'
}

$usedStyles = @()
$unrecognized = @()
$activeSuffixes = @()
$suffixLines = @()

foreach ($lline in $logical) {
    $tokens = $lline -split '\s+'
    $cmd = $tokens[0]

    if ($catMap.ContainsKey($cmd)) {
        $cat = $catMap[$cmd]
        if ($tokens.Count -ge 2) {
            $first = $tokens[1]
            if ($hybridStyles -contains $first) {
                for ($k = 2; $k -lt $tokens.Count; $k++) {
                    $tok = $tokens[$k]
                    $base = (Get-BaseAndSuffix $tok).Base
                    if ($Index[$cat].ContainsKey($base)) {
                        $usedStyles += [pscustomobject]@{ Category = $cat; Style = $tok; Line = $lline }
                    }
                }
            }
            else {
                $base = (Get-BaseAndSuffix $first).Base
                if ($Index[$cat].ContainsKey($base)) {
                    $usedStyles += [pscustomobject]@{ Category = $cat; Style = $first; Line = $lline }
                }
                else {
                    $unrecognized += [pscustomobject]@{ Category = $cat; Style = $first; Line = $lline }
                }
            }
        }
    }
    elseif ($cmd -eq 'fix' -and $tokens.Count -ge 4) {
        $base = (Get-BaseAndSuffix $tokens[3]).Base
        if ($Index['Fix'].ContainsKey($base)) {
            $usedStyles += [pscustomobject]@{ Category = 'Fix'; Style = $tokens[3]; Line = $lline }
        }
        else {
            $unrecognized += [pscustomobject]@{ Category = 'Fix'; Style = $tokens[3]; Line = $lline }
        }
    }
    elseif ($cmd -eq 'compute' -and $tokens.Count -ge 4) {
        $base = (Get-BaseAndSuffix $tokens[3]).Base
        if ($Index['Compute'].ContainsKey($base)) {
            $usedStyles += [pscustomobject]@{ Category = 'Compute'; Style = $tokens[3]; Line = $lline }
        }
        else {
            $unrecognized += [pscustomobject]@{ Category = 'Compute'; Style = $tokens[3]; Line = $lline }
        }
    }
    elseif ($cmd -eq 'suffix' -and $tokens.Count -ge 2) {
        $arg = $tokens[1]
        switch ($arg) {
            'off' { $activeSuffixes = @() }
            'on' { }
            'hybrid' { $activeSuffixes = @($tokens[2..($tokens.Count - 1)]) }
            default { $activeSuffixes = @($arg) }
        }
        $suffixLines += $lline
    }
}

# ---------------------------------------------------------------------------
# 4. Cross-check which accelerator packages are enabled in the chosen build
# ---------------------------------------------------------------------------
$buildStatus = [ordered]@{}
$cmakeCache = Join-Path $RepoRoot "build\$BuildDir\CMakeCache.txt"
if (Test-Path $cmakeCache) {
    foreach ($line in Get-Content $cmakeCache) {
        if ($line -match '^PKG_([A-Z]+):BOOL=(ON|OFF)$') {
            $pkgName = $Matches[1]
            $val = $Matches[2]
            foreach ($v in $VariantKeys) {
                if ($PkgNames[$v] -eq $pkgName) { $buildStatus[$v] = $val }
            }
        }
    }
}

# ---------------------------------------------------------------------------
# 5. Build per-style report items
# ---------------------------------------------------------------------------
$reportItems = @()
foreach ($used in $usedStyles) {
    $bs = Get-BaseAndSuffix $used.Style
    $baseName = $bs.Base
    $directVariant = $bs.Suffix
    $entry = $Index[$used.Category][$baseName]

    $available = @($VariantKeys | Where-Object { $entry[$_] })

    $effective = @()
    if ($directVariant) { $effective += $directVariant }
    $effective += $activeSuffixes

    $activeForThis = @($available | Where-Object { $effective -contains $_ })
    $unusedAvailable = @($available | Where-Object { $effective -notcontains $_ })

    $reportItems += [pscustomobject]@{
        Category       = $used.Category
        BaseName       = $baseName
        UsedAs         = $used.Style
        Line           = $used.Line
        Available      = $available
        ActiveVariants = $activeForThis
        UnusedVariants = $unusedAvailable
    }
}

$totalWithAvail = @($reportItems | Where-Object { $_.Available.Count -gt 0 }).Count
$coverage = @{}
foreach ($v in $VariantKeys) {
    $coverage[$v] = @($reportItems | Where-Object { $_.UnusedVariants -contains $v }).Count
}
$maxCov = 0
if ($coverage.Values.Count -gt 0) { $maxCov = ($coverage.Values | Measure-Object -Maximum).Maximum }
$bestSuffixes = @($VariantKeys | Where-Object { $coverage[$_] -eq $maxCov -and $maxCov -gt 0 })

$commit = $null
if (Test-Path (Join-Path $RepoRoot "lammps-src\.git")) {
    $commit = (& git -C (Join-Path $RepoRoot "lammps-src") rev-parse --short HEAD 2>$null)
}

# ---------------------------------------------------------------------------
# 6. Print report
# ---------------------------------------------------------------------------
Write-Host ""
Write-Host "=== Accelerator coverage report: $InputPath ==="
if ($commit) { Write-Host "LAMMPS source commit: $commit" }

if ($activeSuffixes.Count -gt 0) {
    Write-Host "Active suffix in script: $($activeSuffixes -join ', ')"
}
else {
    Write-Host "Active suffix in script: (none)"
}

if ($buildStatus.Count -gt 0) {
    $bsStr = ($VariantKeys | ForEach-Object {
            $v = if ($buildStatus.Contains($_)) { $buildStatus[$_] } else { '?' }
            "$($PkgNames[$_])=$v"
        }) -join ', '
    Write-Host "Build status (build/$BuildDir): $bsStr"
}
else {
    Write-Host "Build status: build/$BuildDir/CMakeCache.txt not found (run Build-Lammps.ps1 first to check)"
}

Write-Host ""
Write-Host "Styles used in this script:"
foreach ($item in $reportItems) {
    Write-Host ""
    Write-Host ("  [{0}] {1}" -f $DisplayLabel[$item.Category], $item.BaseName)
    Write-Host ("    {0}" -f $item.Line)
    if ($item.Available.Count -eq 0) {
        Write-Host "    No accelerated variant available for this style."
        continue
    }
    Write-Host ("    Accelerated variants available: {0}" -f ($item.Available -join ', '))
    if ($item.ActiveVariants.Count -gt 0) {
        Write-Host ("    -> already accelerated via: {0}" -f ($item.ActiveVariants -join ', '))
    }
    if ($item.UnusedVariants.Count -gt 0) {
        Write-Host ("    -> NOT currently using: {0}" -f ($item.UnusedVariants -join ', '))
        foreach ($v in $item.UnusedVariants) {
            $enabled = if ($buildStatus.Contains($v)) { $buildStatus[$v] } else { 'unknown' }
            Write-Host ("       suggest: -sf $v   (needs PKG_$($PkgNames[$v]); build/$BuildDir has it $enabled)")
        }
    }
}

if ($unrecognized.Count -gt 0) {
    Write-Host ""
    Write-Host "Commands not matched to a known style (skipped):"
    foreach ($u in $unrecognized) {
        Write-Host ("  [{0}] {1} -- {2}" -f $DisplayLabel[$u.Category], $u.Style, $u.Line)
    }
}

Write-Host ""
Write-Host "Summary:"
Write-Host ("  {0} style command(s) recognized, {1} have at least one accelerated variant." -f $reportItems.Count, $totalWithAvail)
if ($bestSuffixes.Count -gt 0) {
    foreach ($v in $bestSuffixes) {
        Write-Host ("  '-sf {0}' (or 'suffix {0}') would newly accelerate {1}/{2} eligible style command(s)." -f $v, $coverage[$v], $totalWithAvail)
        Write-Host ("     run with: $($ExtraArgsHint[$v])")
        $enabled = if ($buildStatus.Contains($v)) { $buildStatus[$v] } else { $null }
        if ($enabled -ne 'ON') {
            Write-Host ("     requires rebuild with -DPKG_$($PkgNames[$v])=ON, e.g.:")
            Write-Host ("       .\tools\bench\Build-Lammps.ps1 -BuildDir $v -Packages $($PkgNames[$v]) -Reconfigure")
        }
    }
}
elseif ($totalWithAvail -gt 0) {
    Write-Host "  All available accelerated variants for the used styles are already active."
}
else {
    Write-Host "  No accelerated variants available for any recognized style in this script."
}
Write-Host ""

# ---------------------------------------------------------------------------
# 7. Optional JSON output
# ---------------------------------------------------------------------------
if ($Json) {
    $jsonObj = [ordered]@{
        InputFile      = $InputPath
        LammpsCommit   = $commit
        BuildDir       = $BuildDir
        BuildStatus    = $buildStatus
        ActiveSuffixes = $activeSuffixes
        Items          = $reportItems
        Unrecognized   = $unrecognized
        Coverage       = $coverage
        BestSuffixes   = $bestSuffixes
    }
    $jsonDir = Split-Path -Parent $Json
    if ($jsonDir -and -not (Test-Path $jsonDir)) { New-Item -ItemType Directory -Force -Path $jsonDir | Out-Null }
    $jsonObj | ConvertTo-Json -Depth 10 | Out-File -FilePath $Json -Encoding utf8
    Write-Host "Saved JSON report: $Json"
}
