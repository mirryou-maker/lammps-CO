<#
.SYNOPSIS
Configure and build LAMMPS from lammps-src/cmake using CMake + MSVC.

.DESCRIPTION
Out-of-source build into <repo>\build\<BuildDir> using the Visual Studio
generator (multi-config). Re-runs `cmake` configure only if the build
directory has no CMakeCache.txt yet, or -Reconfigure is passed.

.EXAMPLE
.\Build-Lammps.ps1 -BuildDir serial -Config Release

.EXAMPLE
.\Build-Lammps.ps1 -BuildDir omp -Config Release -Packages OPENMP -Reconfigure

.EXAMPLE
.\Build-Lammps.ps1 -BuildDir serial-avx2 -Config Release -ExtraFlags "/arch:AVX2" -Reconfigure

.EXAMPLE
.\Build-Lammps.ps1 -BuildDir serial-lto -Config Release -Ipo -Reconfigure
#>

[CmdletBinding()]
param(
    [string]$BuildDir = "serial",
    [string]$Config = "Release",
    [string]$Generator = "Visual Studio 17 2022",
    [string[]]$Packages = @(),
    [int]$Jobs = 0,
    [string]$ExtraFlags = "",
    [switch]$Ipo,
    [switch]$Reconfigure,
    # Path to the LAMMPS source tree (default: lammps-src/ next to repo root).
    # Override if you cloned LAMMPS under a different name:
    #   .\Build-Lammps.ps1 -LammpsDir C:\path\to\lammps
    [string]$LammpsDir = ""
)

$ErrorActionPreference = "Stop"

$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path

# Resolve LAMMPS source directory
if ($LammpsDir -eq "") {
    $LammpsDir = Join-Path $RepoRoot "lammps-src"
}
if (-not (Test-Path $LammpsDir)) {
    Write-Error ("LAMMPS source directory not found: '$LammpsDir'`n" +
        "Clone LAMMPS as 'lammps-src' beside this repo, or pass -LammpsDir <path>.`n" +
        "  git clone https://github.com/lammps/lammps.git '$LammpsDir'")
    exit 1
}
$SourceDir = Join-Path $LammpsDir "cmake"
$BuildPath = Join-Path $RepoRoot "build\$BuildDir"

if ($Jobs -le 0) { $Jobs = [Environment]::ProcessorCount }

New-Item -ItemType Directory -Force -Path $BuildPath | Out-Null

$cacheFile = Join-Path $BuildPath "CMakeCache.txt"
if ($Reconfigure -or -not (Test-Path $cacheFile)) {
    $cmakeArgs = @("-S", $SourceDir, "-B", $BuildPath, "-G", $Generator, "-A", "x64")
    foreach ($pkg in $Packages) {
        $cmakeArgs += "-DPKG_$pkg=ON"
    }
    if ($ExtraFlags) {
        $cmakeArgs += "-DCMAKE_CXX_FLAGS=$ExtraFlags"
        $cmakeArgs += "-DCMAKE_C_FLAGS=$ExtraFlags"
    }
    if ($Ipo) {
        $cmakeArgs += "-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON"
    }
    Write-Host "Configuring: cmake $($cmakeArgs -join ' ')"
    & cmake @cmakeArgs
    if ($LASTEXITCODE -ne 0) { throw "cmake configure failed with exit code $LASTEXITCODE" }
}
else {
    Write-Host "Using existing configuration in $BuildPath (pass -Reconfigure to redo)"
}

Write-Host "Building ($Config, --parallel $Jobs) ..."
& cmake --build $BuildPath --config $Config --parallel $Jobs
if ($LASTEXITCODE -ne 0) { throw "cmake build failed with exit code $LASTEXITCODE" }

$exe = Join-Path $BuildPath "$Config\lmp.exe"
if (Test-Path $exe) {
    Write-Host "Build OK: $exe"
}
else {
    Write-Warning "Build finished but $exe was not found"
}
