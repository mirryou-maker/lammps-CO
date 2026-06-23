#!/usr/bin/env bash
# apply_patch.sh — Apply LLM-optimized LAMMPS files to a LAMMPS source tree
#
# Usage:  bash apply_patch.sh <path-to-lammps-root>
# Example: bash apply_patch.sh ~/software/lammps
#
# What it does:
#   1. Copies 40 new/updated pair_*_omp.{h,cpp} into lammps/src/OPENMP/
#   2. Copies 15 A-3 optimized pair_*.cpp into lammps/src/
#   3. Copies 4 A-4 nm/cut pair_*.cpp into lammps/src/EXTRA-PAIR/
#   4. No CMakeLists.txt changes needed — LAMMPS auto-discovers _omp.cpp files.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

LAMMPS_ROOT="${1:-}"
if [[ -z "$LAMMPS_ROOT" ]]; then
  echo "Usage: $0 <path-to-lammps-root>"
  echo "  e.g.: $0 ~/software/lammps"
  exit 1
fi

if [[ ! -d "$LAMMPS_ROOT/src" ]]; then
  echo "ERROR: $LAMMPS_ROOT/src not found. Is this a LAMMPS source tree?"
  exit 1
fi

echo "=== Applying LLM-optimized LAMMPS patches to: $LAMMPS_ROOT ==="
echo ""

# ── A-1: Copy new OMP pair_style files ──────────────────────────────────────
OMP_SRC="$REPO_ROOT/src/OPENMP"
OMP_DST="$LAMMPS_ROOT/src/OPENMP"

if [[ ! -d "$OMP_DST" ]]; then
  echo "WARNING: $OMP_DST not found. Is PKG_OPENMP available in this LAMMPS?"
  echo "         Skipping OMP files."
else
  count=0
  for f in "$OMP_SRC"/pair_*_omp.cpp "$OMP_SRC"/pair_*_omp.h; do
    fname=$(basename "$f")
    if [[ -e "$OMP_DST/$fname" ]]; then
      echo "  [SKIP] $fname already exists in OPENMP/"
    else
      cp "$f" "$OMP_DST/$fname"
      echo "  [COPY] OPENMP/$fname"
      ((count++))
    fi
  done
  echo "  → $count new OMP files copied."
fi

echo ""

# ── A-3: Copy optimized standard pair files ─────────────────────────────────
STD_SRC="$REPO_ROOT/src"
STD_DST="$LAMMPS_ROOT/src"

a3_files=(
  pair_born.cpp pair_buck.cpp pair_buck_coul_cut.cpp
  pair_coul_cut.cpp pair_coul_debye.cpp pair_coul_dsf.cpp pair_coul_wolf.cpp
  pair_lj_cut.cpp pair_lj_cut_coul_cut.cpp pair_lj_expand.cpp
  pair_morse.cpp pair_soft.cpp pair_table.cpp pair_yukawa.cpp pair_zbl.cpp
)

count=0
for f in "${a3_files[@]}"; do
  if [[ ! -e "$STD_SRC/$f" ]]; then
    echo "  [WARN] $f not found in repo src/ — skipping"
    continue
  fi
  cp "$STD_SRC/$f" "$STD_DST/$f"
  echo "  [COPY] src/$f (A-3 optimized)"
  ((count++))
done
echo "  → $count standard pair files replaced."

echo ""

# ── A-4: Copy nm/cut optimized EXTRA-PAIR files ──────────────────────────────
EP_SRC="$REPO_ROOT/src/EXTRA-PAIR"
EP_DST="$LAMMPS_ROOT/src/EXTRA-PAIR"

a4_files=(
  pair_nm_cut.cpp pair_nm_cut_coul_cut.cpp
  pair_nm_cut_coul_long.cpp pair_nm_cut_split.cpp
)

if [[ ! -d "$EP_DST" ]]; then
  echo "WARNING: $EP_DST not found. Build with PKG_EXTRA-PAIR=yes."
  echo "         Skipping A-4 nm/cut files."
else
  count=0
  for f in "${a4_files[@]}"; do
    if [[ ! -e "$EP_SRC/$f" ]]; then
      echo "  [WARN] $f not found in repo src/EXTRA-PAIR/ — skipping"
      continue
    fi
    cp "$EP_SRC/$f" "$EP_DST/$f"
    echo "  [COPY] EXTRA-PAIR/$f (A-4 pow reduction)"
    ((count++))
  done
  echo "  → $count nm/cut files replaced."
fi

echo ""
echo "=== Patch applied successfully ==="
echo ""
echo "Next step: build LAMMPS (see INSTALL.md)"
