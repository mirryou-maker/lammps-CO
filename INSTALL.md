# Installation Guide — LLM-Optimized LAMMPS

This guide explains how to download, patch, and build the LLM-optimized LAMMPS
(OpenMP A-1 variants + A-3 restrict/fxtmp backport + A-4 nm/cut pow() reduction)
on Linux/macOS and Windows.

---

## Overview of changes

| Category | Files | Effect |
|---|---|---|
| **A-1 — 36 new OMP variants** | `src/OPENMP/pair_*_omp.{h,cpp}` (80 files) | Adds `pair_style xxx/omp` for 36 previously unparallelized styles → 4–5× speedup on multi-core nodes |
| **A-3 — restrict backport** | `src/pair_*.cpp` (15 files) | Enables compiler auto-vectorization of standard (serial) hotloops → +4–6% serial speed |
| **A-4 — pow() reduction (nm/cut)** | `src/EXTRA-PAIR/pair_nm_cut*.cpp` (4 files) + `src/OPENMP/pair_nm_cut*_omp.cpp` (4 files) | Replaces 4 `pow()` calls per pair with 2 in force block → **−39.8% loop time** for nm/cut family |
| **A-4 — pow() reduction (lj/pirani)** | `src/EXTRA-PAIR/pair_lj_pirani.cpp`, `src/OPENMP/pair_lj_pirani_omp.cpp` | Eliminates 2 redundant `pow()` in energy block → **−23.4% loop time** when energy computed every step |
| **A-4 — pow() reduction (mie/cut)** | `src/EXTRA-PAIR/pair_mie_cut.cpp` | Eliminates redundant `pow()` in compute_outer (RESPA multi-timestep) |

**No CMakeLists.txt changes required.** LAMMPS auto-discovers `_omp.cpp` files via `RegisterStylesExt`.

---

## Step 1 — Clone LAMMPS (develop branch)

```bash
git clone --depth=1 --branch=develop https://github.com/lammps/lammps.git
```

> **Minimum version**: commit `91d4111` or later (2025-06). The patches were
> developed against this commit.

---

## Step 2 — Download this repository

```bash
# Option A: git clone (requires GitHub account with repo access)
git clone https://github.com/mirryou-maker/lammps-llm-omp-optimization.git

# Option B: Download ZIP from GitHub Release page
#   https://github.com/mirryou-maker/lammps-llm-omp-optimization/releases/latest
#   → "Source code (zip)" → unzip
```

---

## Step 3 — Apply the patch

### Linux / macOS

```bash
cd lammps-llm-omp-optimization
bash scripts/apply_patch.sh /path/to/lammps
```

### Windows (PowerShell)

```powershell
cd lammps-llm-omp-optimization
.\scripts\apply_patch.ps1 -LammpsRoot C:\path\to\lammps
```

The script:
- Copies 80 new/updated OMP files into `lammps/src/OPENMP/`
- Replaces 15 standard pair files in `lammps/src/` with A-3 optimized versions
- Copies 4 A-4 nm/cut files into `lammps/src/EXTRA-PAIR/`
- Skips files that already exist (safe to re-run)

### Manual copy (if scripts are unavailable)

```bash
# Copy new OMP pair files
cp src/OPENMP/*.cpp  /path/to/lammps/src/OPENMP/
cp src/OPENMP/*.h    /path/to/lammps/src/OPENMP/

# Copy A-3 optimized standard files
cp src/pair_*.cpp    /path/to/lammps/src/

# Copy A-4 nm/cut optimized files (requires PKG_EXTRA-PAIR)
cp src/EXTRA-PAIR/*.cpp  /path/to/lammps/src/EXTRA-PAIR/
```

---

## Step 4 — Build LAMMPS

### Linux / macOS (CMake — recommended)

```bash
cd /path/to/lammps
mkdir build && cd build

cmake ../cmake \
  -D CMAKE_BUILD_TYPE=Release \
  -D BUILD_SHARED_LIBS=no \
  -D PKG_OPENMP=yes \
  -D PKG_KSPACE=yes \
  -D PKG_MOLECULE=yes \
  -D PKG_EXTRA-PAIR=yes \
  -D CMAKE_CXX_FLAGS="-O3 -march=native -funroll-loops"

make -j$(nproc)
```

> **Tip — FEP/alchemical simulations** also require `PKG_FEP`:
> ```bash
> cmake ../cmake ... -D PKG_FEP=yes -D PKG_OPENMP=yes
> ```

> **Tip — AVX2 + fast-math** (best performance, non-IEEE FP):
> ```bash
> cmake ../cmake ... -D CMAKE_CXX_FLAGS="-O3 -march=native -ffp-contract=fast"
> ```

### Windows (CMake + MSVC)

```powershell
cd C:\path\to\lammps
mkdir build; cd build

cmake ..\cmake `
  -G "Visual Studio 17 2022" -A x64 `
  -D CMAKE_BUILD_TYPE=Release `
  -D PKG_OPENMP=yes `
  -D PKG_KSPACE=yes `
  -D PKG_MOLECULE=yes `
  -D PKG_EXTRA-PAIR=yes

cmake --build . --config Release --parallel
```

### HPC clusters (Makefile-based build)

```bash
cd /path/to/lammps/src

# Install the OPENMP package and required packages
make yes-OPENMP yes-KSPACE yes-MOLECULE yes-EXTRA-PAIR

# Copy the new OMP files (they are already in src/OPENMP/)
# Build with your cluster's MPI + OpenMP compiler wrapper
make -j8 mpi LMP_INC="-fopenmp" CC=mpicxx CCFLAGS="-O3 -march=native"
```

---

## Step 5 — Verify the build

```bash
# Check the binary exists and list OMP pair styles
./lmp -h 2>&1 | grep "/omp" | grep "^pair"
```

You should see the 36 new styles listed, e.g.:
```
pair_style born/coul/dsf/omp
pair_style born/coul/dsf/cs/omp
pair_style coul/ctip/omp
pair_style lj/charmmfsw/coul/charmmfsh/omp
pair_style lj/class2/soft/omp
pair_style morse/soft/omp
pair_style nm/cut/split/omp
...
```

Run the bundled melt benchmark to confirm speedup:

```bash
# Serial baseline
./lmp -in ../examples/melt/in.melt -log serial.log

# OMP 4 threads
OMP_NUM_THREADS=4 ./lmp -sf omp -in ../examples/melt/in.melt -log omp4.log

# Compare "Loop time" in the two logs — expect ~4× improvement
grep "Loop time" serial.log omp4.log
```

---

## Step 6 — Using the new pair styles

### Drop-in replacement (suffix)

The easiest way to use the new OMP variants is the `-sf omp` command-line flag,
which automatically appends `/omp` to all pair, bond, angle, etc. styles:

```bash
OMP_NUM_THREADS=8 ./lmp -sf omp -pk omp 8 -in in.your_simulation
```

This works for **all** 36 newly added styles with no script changes.

### Explicit style selection

```lammps
# in LAMMPS input script
package omp 4

# Use the OMP variant directly:
pair_style  lj/charmmfsw/coul/charmmfsh/omp  12.0 14.0
pair_style  coul/ctip/omp  10.0
pair_style  morse/soft/omp  2  0.5  12.0

# For FEP alchemical calculations:
pair_style  lj/cut/coul/long/soft/omp  2  0.5  12.0
pair_coeff  * * ...
fix         fep all adapt/fep ...
```

### A-3 optimized serial styles (no script changes needed)

The 15 modified `pair_*.cpp` files are drop-in replacements for the originals.
They produce **bit-identical results** and require no changes to LAMMPS input scripts.
Benefits are automatic when using the standard (non-OMP) styles:

```lammps
pair_style  lj/cut  2.5       # already uses A-3 optimized code after patching
pair_style  lj/cut/coul/cut  10.0
pair_style  morse  10.0
```

---

## New pair styles added (36 total)

| Style name (`/omp` suffix) | Base package |
|---|---|
| `born/coul/dsf/omp` | KSPACE |
| `born/coul/dsf/cs/omp` | KSPACE |
| `born/coul/long/cs/omp` | KSPACE |
| `born/coul/wolf/cs/omp` | KSPACE |
| `bpm/spring/omp` | BPM |
| `buck/coul/long/cs/omp` | KSPACE |
| `coul/ctip/omp` | (base) |
| `coul/cut/dielectric/omp` | DIELECTRIC |
| `coul/cut/soft/gapsys/omp` | FEP |
| `coul/esp/omp` | (base) |
| `coul/exclude/omp` | (base) |
| `coul/long/cs/omp` | KSPACE |
| `coul/long/dielectric/omp` | DIELECTRIC |
| `coul/slater/cut/omp` | (base) |
| `coul/slater/long/omp` | KSPACE |
| `coul/tt/omp` | (base) |
| `coul/wolf/cs/omp` | KSPACE |
| `dpd/fdt/omp` | DPD-MESO |
| `kolmogorov/crespi/z/omp` | INTERLAYER |
| `lebedeva/z/omp` | INTERLAYER |
| `lj/charmmfsw/coul/charmmfsh/omp` | KSPACE |
| `lj/charmmfsw/coul/long/omp` | KSPACE |
| `lj/class2/coul/cut/soft/omp` | FEP |
| `lj/class2/coul/long/cs/omp` | KSPACE |
| `lj/class2/coul/long/soft/omp` | FEP |
| `lj/class2/soft/omp` | FEP |
| `lj/cut/coul/esp/omp` | (base) |
| `lj/cut/coul/long/cs/omp` | KSPACE |
| `lj/cut/coul/msm/dielectric/omp` | DIELECTRIC |
| `lj/cut/dipole/long/omp` | DIPOLE |
| `lj/expand/coul/long/omp` | KSPACE |
| `mdpd/omp` | DPD-MESO |
| `morse/soft/omp` | FEP |
| `nm/cut/split/omp` | (base) |
| `rheo/solid/omp` | RHEO |
| `thole/omp` | DRUDE |

---

## A-3 optimized standard pair files

The following 15 files have been optimized with `__restrict__` / `_noalias` hints
and `fxtmp/fytmp/fztmp` register-accumulator pattern:

`pair_born`, `pair_buck`, `pair_buck_coul_cut`, `pair_coul_cut`, `pair_coul_debye`,
`pair_coul_dsf`, `pair_coul_wolf`, `pair_lj_cut`, `pair_lj_cut_coul_cut`,
`pair_lj_expand`, `pair_morse`, `pair_soft`, `pair_table`, `pair_yukawa`, `pair_zbl`

**Performance**: +4.2% serial loop-time improvement on `pair_lj_cut` (32,000 atoms);
results are **bit-identical** to the unoptimized originals.

---

## A-4 pow() reduction (nm/cut, lj/pirani, mie/cut — EXTRA-PAIR)

The following files replace original EXTRA-PAIR/OPENMP implementations to eliminate
redundant transcendental `pow()` calls per neighbor pair:

### nm/cut family (4→2 pow per pair, force block — every step)

| File | Optimization |
|---|---|
| `src/EXTRA-PAIR/pair_nm_cut.cpp` | A-3 + A-4 |
| `src/EXTRA-PAIR/pair_nm_cut_coul_cut.cpp` | A-3 + A-4 |
| `src/EXTRA-PAIR/pair_nm_cut_coul_long.cpp` | A-3 + A-4 |
| `src/EXTRA-PAIR/pair_nm_cut_split.cpp` | A-3 + A-4 + eflag dedup |

**Performance**: **−39.8–40.9% loop time** (1.65–1.69×), 864–6,912 atoms, serial 1T.

### lj/pirani (2 pow eliminated per pair, energy block — when eflag active)

| File | Optimization |
|---|---|
| `src/EXTRA-PAIR/pair_lj_pirani.cpp` | A-4 eflag block pow reuse |
| `src/OPENMP/pair_lj_pirani_omp.cpp` | A-4 EFLAG template block pow reuse |

**Performance**: **−23.4% loop time** (1.31×) when energy is computed every step
(`thermo 1` or energy minimization). Benefit scales with thermo frequency.

### mie/cut (compute_outer RESPA fix)

| File | Optimization |
|---|---|
| `src/EXTRA-PAIR/pair_mie_cut.cpp` | A-4 compute_outer: lift r2inv/rgamA/rgamR before conditional blocks |

Eliminates up to 4 redundant `pow()` calls per pair in the RESPA `compute_outer()`
path (multi-timestep integration). Standard `compute()` was already optimal.

**Requires**: `PKG_EXTRA-PAIR=yes` in the cmake build for all EXTRA-PAIR files.

---

## Performance expectations

| Configuration | Expected speedup vs. serial |
|---|---|
| OMP 4 threads (`-sf omp -pk omp 4`) | ~4× (Pair time; Amdahl-limited) |
| OMP 8 threads (`-sf omp -pk omp 8`) | ~5× |
| A-3 only (no OMP) | ~4–6% |
| A-4 nm/cut (pow reduction) | **−39.8% loop time** (1.65×, any thread count) |
| OMP 8t + AVX2 + A-3 | ~5–6× |

> Tested on: AMD Ryzen 9 / 32,000-atom LJ melt benchmark (OMP/A-3) and 864–6,912
> atom nm/cut benchmark (A-4), LAMMPS develop branch commit `91d4111`, MSVC 2022.
> Actual performance varies by hardware and system size.

---

## Troubleshooting

**Q: Build fails with "unknown pair style xxx/omp"**
→ The style requires a specific package (`PKG_KSPACE`, `PKG_FEP`, etc.). See the
"Base package" column in the table above.

**Q: `lmp -h` doesn't show new styles after patching**
→ Confirm the `.cpp` and `.h` files are in `lammps/src/OPENMP/` and rebuild from scratch:
```bash
rm -rf build && mkdir build && cd build && cmake ... && make -j$(nproc)
```

**Q: Numerical results differ from unpatched LAMMPS**
→ The A-3 files are bit-identical for IEEE-754 arithmetic. If using `-ffp-contract=fast`
or `-ffast-math`, minor floating-point reordering differences (< 1 ULP) are expected
and do not affect physics.

**Q: Windows MSVC build: OpenMP not found**
→ Install "C++ OpenMP" via Visual Studio Installer → modify → Individual components.
Then add `-D OpenMP_CXX_FLAGS="/openmp"` to the cmake command.

---

## Citation

If you use this optimized LAMMPS in your research, please cite:

> You, C.-Y. "AI-Assisted Systematic Optimization of the LAMMPS Molecular Dynamics
> Simulator via Large Language Model Coding Agents."
> *npj Computational Materials* (2026). https://github.com/mirryou-maker/lammps-llm-omp-optimization

Also cite the original LAMMPS paper:

> Thompson, A. P. et al. LAMMPS — A flexible simulation tool for particle-based
> materials modeling. *Comput. Phys. Commun.* **271**, 108171 (2022).
