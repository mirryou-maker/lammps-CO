#!/usr/bin/env bash
# apply_patch.sh — Apply LLM-optimized LAMMPS files to a LAMMPS source tree
#
# Usage:  bash apply_patch.sh <path-to-lammps-root>
# Example: bash apply_patch.sh ~/software/lammps
#
# What it does:
#   1. Copies new/updated pair_*_omp.{h,cpp} into lammps/src/OPENMP/        (A-1)
#   2. Copies every A-3 optimized pair_*.cpp into its LAMMPS package dir    (A-3)
#      (src/ plus 37 package subdirectories — inventory in
#       Paper/tableS6_a3_files.csv)
#   3. A-4 pow()-reduction files ship inside the same tree (src/EXTRA-PAIR/,
#      src/OPENMP/) and are copied by the same walk.
#   4. No CMakeLists.txt changes needed — LAMMPS auto-discovers _omp.cpp files.
#
# A-3 files are copied by walking this repository's src/ tree instead of from a
# hard-coded list, so the script cannot drift out of sync with the shipped code.
# Only files that already exist in the target tree are replaced, so a package
# the user did not enable is reported rather than silently created.

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

# ── A-1: new / updated OMP pair styles ──────────────────────────────────────
OMP_SRC="$REPO_ROOT/src/OPENMP"
OMP_DST="$LAMMPS_ROOT/src/OPENMP"

echo "-- A-1: OpenMP pair styles --"
if [[ ! -d "$OMP_DST" ]]; then
  echo "WARNING: $OMP_DST not found. Is PKG_OPENMP available in this LAMMPS?"
  echo "         Skipping OMP files."
else
  new=0; upd=0
  for f in "$OMP_SRC"/pair_*_omp.cpp "$OMP_SRC"/pair_*_omp.h; do
    [[ -e "$f" ]] || continue
    fname=$(basename "$f")
    if [[ -e "$OMP_DST/$fname" ]]; then
      upd=$((upd + 1))
    else
      new=$((new + 1))
    fi
    cp "$f" "$OMP_DST/$fname"
  done
  echo "  -> $new new, $upd updated OMP files."
fi

echo ""

# ── A-3 / A-4: restrict + local-accumulator (and pow-reduced) pair files ────
# Layout mirrors LAMMPS:
#   src/pair_x.cpp           -> <lammps>/src/pair_x.cpp
#   src/<PACKAGE>/pair_x.cpp -> <lammps>/src/<PACKAGE>/pair_x.cpp
echo "-- A-3/A-4: restrict/accumulator + pow()-reduced pair files --"
copied=0; missing=0
skipped_pkg=""

while IFS= read -r f; do
  rel="${f#"$REPO_ROOT/src/"}"          # MANYBODY/pair_eam.cpp | pair_lj_cut.cpp
  case "$rel" in
    OPENMP/*) continue ;;               # handled above
  esac
  dst="$LAMMPS_ROOT/src/$rel"
  dstdir="$(dirname "$dst")"
  if [[ ! -d "$dstdir" ]]; then
    skipped_pkg="$skipped_pkg $(dirname "$rel")"
    continue
  fi
  if [[ ! -e "$dst" ]]; then
    echo "  [WARN] $rel not present in target tree - skipping"
    missing=$((missing + 1))
    continue
  fi
  cp "$f" "$dst"
  copied=$((copied + 1))
done < <(find "$REPO_ROOT/src" -name 'pair_*.cpp' | sort)

echo "  -> $copied pair files replaced."
[[ $missing -gt 0 ]] && echo "  -> $missing not found in target tree."
if [[ -n "$skipped_pkg" ]]; then
  uniq_pkg=$(printf '%s\n' $skipped_pkg | sort -u | tr '\n' ' ')
  echo "  -> packages absent from target tree (not enabled?): $uniq_pkg"
fi

echo ""
echo "=== Patch applied successfully ==="
echo ""
echo "Next step: build LAMMPS (see INSTALL.md), e.g."
echo "  cd $LAMMPS_ROOT && mkdir -p build && cd build"
echo "  cmake ../cmake -D CMAKE_BUILD_TYPE=Release -D PKG_OPENMP=on \\"
echo "        -D BUILD_OMP=on -D BUILD_MPI=on"
echo "  cmake --build . -j"
