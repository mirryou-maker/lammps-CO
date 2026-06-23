# apply_patch.ps1 — Apply LLM-optimized LAMMPS files to a LAMMPS source tree (Windows)
#
# Usage:  .\apply_patch.ps1 -LammpsRoot <path-to-lammps-root>
# Example: .\apply_patch.ps1 -LammpsRoot C:\software\lammps
#
# What it does:
#   1. Copies 40 new/updated pair_*_omp.{h,cpp} into lammps\src\OPENMP\
#   2. Copies 15 A-3 optimized pair_*.cpp into lammps\src\
#   3. Copies 4 A-4 nm/cut pair_*.cpp into lammps\src\EXTRA-PAIR\
#   No CMakeLists.txt changes needed — LAMMPS auto-discovers _omp.cpp files.

param(
  [Parameter(Mandatory=$true)]
  [string]$LammpsRoot
)

$ScriptDir  = Split-Path -Parent $MyInvocation.MyCommand.Path
$RepoRoot   = Split-Path -Parent $ScriptDir

if (-not (Test-Path "$LammpsRoot\src")) {
  Write-Error "ERROR: $LammpsRoot\src not found. Is this a LAMMPS source tree?"
  exit 1
}

Write-Host "=== Applying LLM-optimized LAMMPS patches to: $LammpsRoot ===" -ForegroundColor Cyan
Write-Host ""

# A-1: New OMP pair_style files
$OmpSrc = "$RepoRoot\src\OPENMP"
$OmpDst = "$LammpsRoot\src\OPENMP"

if (-not (Test-Path $OmpDst)) {
  Write-Warning "  $OmpDst not found — skipping OMP files."
} else {
  $count = 0
  Get-ChildItem "$OmpSrc\pair_*_omp.cpp", "$OmpSrc\pair_*_omp.h" | ForEach-Object {
    $dst = "$OmpDst\$($_.Name)"
    if (Test-Path $dst) {
      Write-Host "  [SKIP] $($_.Name) already exists" -ForegroundColor Yellow
    } else {
      Copy-Item $_.FullName $dst
      Write-Host "  [COPY] OPENMP\$($_.Name)" -ForegroundColor Green
      $count++
    }
  }
  Write-Host "  -> $count new OMP files copied."
}

Write-Host ""

# A-3: Optimized standard pair files
$StdSrc = "$RepoRoot\src"
$StdDst = "$LammpsRoot\src"

$a3files = @(
  'pair_born.cpp','pair_buck.cpp','pair_buck_coul_cut.cpp',
  'pair_coul_cut.cpp','pair_coul_debye.cpp','pair_coul_dsf.cpp','pair_coul_wolf.cpp',
  'pair_lj_cut.cpp','pair_lj_cut_coul_cut.cpp','pair_lj_expand.cpp',
  'pair_morse.cpp','pair_soft.cpp','pair_table.cpp','pair_yukawa.cpp','pair_zbl.cpp'
)

$count = 0
foreach ($f in $a3files) {
  $src = "$StdSrc\$f"
  if (-not (Test-Path $src)) {
    Write-Warning "  $f not found in repo src\ — skipping"
    continue
  }
  Copy-Item $src "$StdDst\$f" -Force
  Write-Host "  [COPY] src\$f  (A-3 optimized)" -ForegroundColor Green
  $count++
}
Write-Host "  -> $count standard pair files replaced."

Write-Host ""

# A-4: nm/cut optimized EXTRA-PAIR files
$EpSrc = "$RepoRoot\src\EXTRA-PAIR"
$EpDst = "$LammpsRoot\src\EXTRA-PAIR"

$a4files = @(
  'pair_nm_cut.cpp','pair_nm_cut_coul_cut.cpp',
  'pair_nm_cut_coul_long.cpp','pair_nm_cut_split.cpp',
  'pair_lj_pirani.cpp','pair_mie_cut.cpp'
)

if (-not (Test-Path $EpDst)) {
  Write-Warning "  $EpDst not found — build with PKG_EXTRA-PAIR=yes. Skipping A-4 files."
} else {
  $count = 0
  foreach ($f in $a4files) {
    $src = "$EpSrc\$f"
    if (-not (Test-Path $src)) {
      Write-Warning "  $f not found in repo src\EXTRA-PAIR\ — skipping"
      continue
    }
    Copy-Item $src "$EpDst\$f" -Force
    Write-Host "  [COPY] EXTRA-PAIR\$f  (A-4 pow reduction)" -ForegroundColor Green
    $count++
  }
  Write-Host "  -> $count EXTRA-PAIR pow-reduction files replaced."
}

Write-Host ""

# A-4 OMP: lj/pirani OMP pow()-reduction (replaces existing OMP file)
$a4ompfiles = @('pair_lj_pirani_omp.cpp')
if (Test-Path $OmpDst) {
  $count = 0
  foreach ($f in $a4ompfiles) {
    $src = "$OmpSrc\$f"
    if (-not (Test-Path $src)) {
      Write-Warning "  $f not found in repo src\OPENMP\ — skipping"
      continue
    }
    Copy-Item $src "$OmpDst\$f" -Force
    Write-Host "  [COPY] OPENMP\$f  (A-4 pow reduction)" -ForegroundColor Green
    $count++
  }
  Write-Host "  -> $count OMP pow-reduction files replaced."
}

Write-Host ""
Write-Host "=== Patch applied successfully ===" -ForegroundColor Cyan
Write-Host ""
Write-Host "Next step: build LAMMPS (see INSTALL.md)"
