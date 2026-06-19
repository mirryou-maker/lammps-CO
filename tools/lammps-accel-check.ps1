<#
.SYNOPSIS
Analyzes a LAMMPS input file and reports unused acceleration opportunities.

.DESCRIPTION
Reads a LAMMPS input file, identifies all styles used (pair_style, bond_style,
angle_style, dihedral_style, improper_style), and checks which LAMMPS acceleration
packages (OPENMP, OPT, INTEL, GPU, KOKKOS) have compiled variants for those styles
that are not being used. Prints a prioritized report with the exact command-line
flags needed to activate each available acceleration.

.PARAMETER InputFile
Path to the LAMMPS input script to analyze.

.PARAMETER LammpsSrc
Path to the lammps src/ directory. Default: ../../lammps-src/src relative to this script.

.EXAMPLE
.\lammps-accel-check.ps1 -InputFile ..\..\lammps-src\examples\melt\in.melt

.EXAMPLE
.\lammps-accel-check.ps1 -InputFile C:\sims\in.polymer -LammpsSrc C:\lammps\src
#>

[CmdletBinding()]
param(
    [Parameter(Mandatory)] [string]$InputFile,
    [string]$LammpsSrc = ""
)

$ErrorActionPreference = "Stop"

$ScriptDir = $PSScriptRoot
if (-not $LammpsSrc) {
    $LammpsSrc = Join-Path $ScriptDir "..\lammps-src\src"
}
$LammpsSrc = (Resolve-Path $LammpsSrc).Path
$InputFile = (Resolve-Path $InputFile).Path

# -----------------------------------------------------------------------
# Package definitions: directory, filename suffix, runtime style suffix,
# recommended command-line flags, priority (lower = faster to get benefit)
# -----------------------------------------------------------------------
$packages = [ordered]@{
    OPENMP = [pscustomobject]@{
        Dir      = Join-Path $LammpsSrc "OPENMP"
        FileSfx  = "_omp"
        StyleSfx = "/omp"
        Flags    = "-sf omp -pk omp <N_THREADS>"
        Priority = 1
        Note     = "shared-memory multithreading; easiest win on any multicore CPU"
    }
    OPT = [pscustomobject]@{
        Dir      = Join-Path $LammpsSrc "OPT"
        FileSfx  = "_opt"
        StyleSfx = "/opt"
        Flags    = "-sf opt"
        Priority = 2
        Note     = "template-based branch elimination; ~5-20% on serial"
    }
    INTEL = [pscustomobject]@{
        Dir      = Join-Path $LammpsSrc "INTEL"
        FileSfx  = "_intel"
        StyleSfx = "/intel"
        Flags    = "-sf intel -pk intel <N_THREADS>"
        Priority = 3
        Note     = "Intel-optimized vectorized kernels; best on Intel CPUs"
    }
    GPU = [pscustomobject]@{
        Dir      = Join-Path $LammpsSrc "GPU"
        FileSfx  = "_gpu"
        StyleSfx = "/gpu"
        Flags    = "-sf gpu -pk gpu <N_GPUS>"
        Priority = 4
        Note     = "CUDA/OpenCL GPU offload; large speedup with GPU hardware"
    }
    KOKKOS = [pscustomobject]@{
        Dir      = Join-Path $LammpsSrc "KOKKOS"
        FileSfx  = "_kokkos"
        StyleSfx = "/kk"
        Flags    = "-k on -sf kk"
        Priority = 5
        Note     = "portable performance abstraction (CPU+GPU); requires Kokkos build"
    }
}

# -----------------------------------------------------------------------
# Parse LAMMPS input file for style declarations
# Returns list of {Category, StyleName, Line} objects
# -----------------------------------------------------------------------
function Get-LammpsStyles([string]$path) {
    $styles = [System.Collections.Generic.List[pscustomobject]]::new()
    $seen   = [System.Collections.Generic.HashSet[string]]::new()

    $categories = @("pair", "bond", "angle", "dihedral", "improper", "kspace")
    $hybridTokens = @("hybrid", "hybrid/overlay", "hybrid/molecular", "hybrid/scaled", "hybrid/weighted")

    $lineNum = 0
    foreach ($raw in (Get-Content $path)) {
        $lineNum++
        $line = ($raw -replace '#.*', '').Trim()
        if (-not $line) { continue }

        foreach ($cat in $categories) {
            $keyword = "${cat}_style"
            if ($line -notmatch "^${keyword}\s+(.+)") { continue }

            $tokens = ($Matches[1].Trim() -split '\s+')
            $first  = $tokens[0].ToLower()

            if ($hybridTokens -contains $first) {
                # hybrid: scan remaining tokens; a new sub-style starts with a letter
                # (numeric tokens are coefficients for the previous sub-style)
                $i = 1
                while ($i -lt $tokens.Length) {
                    if ($tokens[$i] -match '^[a-zA-Z]') {
                        $name = $tokens[$i]
                        $uid  = "${cat}:${name}"
                        if ($seen.Add($uid)) {
                            $styles.Add([pscustomobject]@{
                                Category = $cat; StyleName = $name
                                IsHybrid = $true; Line = $lineNum
                            })
                        }
                    }
                    $i++
                }
            } else {
                $name = $tokens[0]
                $uid  = "${cat}:${name}"
                if ($seen.Add($uid)) {
                    $styles.Add([pscustomobject]@{
                        Category = $cat; StyleName = $name
                        IsHybrid = $false; Line = $lineNum
                    })
                }
            }
            break
        }
    }
    return $styles
}

# -----------------------------------------------------------------------
# For a given {category, styleName}, return which packages have a variant
# -----------------------------------------------------------------------
function Get-AvailablePackages([string]$category, [string]$styleName) {
    # Convert style name to C++ file base:
    #   "lj/cut"          -> "pair_lj_cut"
    #   "harmonic"        -> "bond_harmonic"
    #   "eam/alloy"       -> "pair_eam_alloy"
    $fileBase = "${category}_$($styleName -replace '/', '_')"

    $available = [System.Collections.Generic.List[string]]::new()
    foreach ($pkgName in $packages.Keys) {
        $pkg     = $packages[$pkgName]
        $variant = Join-Path $pkg.Dir "${fileBase}$($pkg.FileSfx).cpp"
        if (Test-Path $variant) {
            $available.Add($pkgName)
        }
    }
    return $available
}

# -----------------------------------------------------------------------
# Determine if the style is already an accelerated variant
# -----------------------------------------------------------------------
function Test-IsAccelerated([string]$styleName) {
    $accelSuffixes = @("/omp", "/opt", "/intel", "/gpu", "/kk")
    foreach ($sfx in $accelSuffixes) {
        if ($styleName.EndsWith($sfx)) { return $true }
    }
    return $false
}

# -----------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------
$styles = @(Get-LammpsStyles $InputFile)

Write-Host ""
Write-Host "LAMMPS Acceleration Opportunity Report"
Write-Host "======================================="
Write-Host "Input : $InputFile"
Write-Host "SrcDir: $LammpsSrc"
Write-Host ""

$opportunities = [System.Collections.Generic.List[pscustomobject]]::new()
$alreadyAccel  = [System.Collections.Generic.List[pscustomobject]]::new()
$noVariant     = [System.Collections.Generic.List[pscustomobject]]::new()

foreach ($s in $styles) {
    if (Test-IsAccelerated $s.StyleName) {
        $alreadyAccel.Add($s)
        continue
    }
    # kspace doesn't follow the same file naming; skip variant check
    if ($s.Category -eq "kspace") {
        $noVariant.Add($s)
        continue
    }
    $avail = @(Get-AvailablePackages $s.Category $s.StyleName)
    if ($avail.Count -gt 0) {
        $opportunities.Add([pscustomobject]@{
            Category  = $s.Category
            StyleName = $s.StyleName
            Line      = $s.Line
            Packages  = $avail
        })
    } else {
        $noVariant.Add($s)
    }
}

# Sort opportunities by highest-priority package available
$opportunities = @($opportunities | Sort-Object {
    ($_.Packages | ForEach-Object { $packages[$_].Priority } | Measure-Object -Minimum).Minimum
})

if ($opportunities.Count -gt 0) {
    Write-Host "[ UNUSED ACCELERATION VARIANTS ]"
    Write-Host "---------------------------------"
    foreach ($o in $opportunities) {
        Write-Host ""
        Write-Host "  $($o.Category)_style $($o.StyleName)  (line $($o.Line))"
        foreach ($pkgName in ($o.Packages | Sort-Object { $packages[$_].Priority })) {
            $pkg = $packages[$pkgName]
            Write-Host "    [$pkgName]  add: $($pkg.Flags)"
            Write-Host "              -> $($pkg.Note)"
        }
    }
} else {
    Write-Host "  All detected styles are already accelerated or have no variants."
}

Write-Host ""
if ($alreadyAccel.Count -gt 0) {
    Write-Host "[ ALREADY ACCELERATED ]"
    foreach ($s in $alreadyAccel) {
        Write-Host "  $($s.Category)_style $($s.StyleName)"
    }
    Write-Host ""
}

if ($noVariant.Count -gt 0) {
    Write-Host "[ NO VARIANT FOUND (standard code only) ]"
    foreach ($s in $noVariant) {
        Write-Host "  $($s.Category)_style $($s.StyleName)"
    }
    Write-Host ""
}

Write-Host "======================================="
Write-Host "Styles detected       : $($styles.Count)"
Write-Host "Already accelerated   : $($alreadyAccel.Count)"
Write-Host "Unused opportunities  : $($opportunities.Count)"
Write-Host "No variant available  : $($noVariant.Count)"
Write-Host ""
if ($opportunities.Count -gt 0) {
    Write-Host "RECOMMENDATION: Add the following to your lammps command line:"
    # Suggest the highest-priority package that covers the most styles
    $pkgCoverage = @{}
    foreach ($o in $opportunities) {
        foreach ($p in $o.Packages) {
            if (-not $pkgCoverage.ContainsKey($p)) { $pkgCoverage[$p] = 0 }
            $pkgCoverage[$p]++
        }
    }
    $best = $pkgCoverage.GetEnumerator() | Sort-Object {
        $packages[$_.Key].Priority * 1000 - $_.Value
    } | Select-Object -First 1
    $bpkg = $packages[$best.Key]
    Write-Host "  $($bpkg.Flags)"
    Write-Host "  (covers $($best.Value) of $($opportunities.Count) style(s) with $($best.Key) package)"
}
