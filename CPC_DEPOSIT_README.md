# LAMMPS-CO

Extended OpenMP coverage and compiler-friendly force kernels for LAMMPS pair styles.

Deposited with the manuscript *"LAMMPS-CO: extended OpenMP coverage and
compiler-friendly force kernels for 261 LAMMPS pair styles"*, submitted to
*Computer Physics Communications* as a Computer Programs in Physics paper.

| | |
|---|---|
| Author | Chun-Yeol You, DGIST (<cyyou@dgist.ac.kr>) |
| License | GPLv2 (matching LAMMPS) |
| Language | C++ (C++11), OpenMP |
| Developed against | LAMMPS `91d4111a9` (develop, 5 June 2026) |
| Repository | <https://github.com/mirryou-maker/lammps-CO> |

## What this is

LAMMPS threads a pair style within an MPI rank only when an OpenMP variant of
that style exists, and fewer than half of its 307 pair styles have one. This
package supplies **36 new OpenMP pair-style variants**, raising coverage to 193
of 307, and rewrites the serial force kernels of **225 further pair files across
38 packages** so the compiler can vectorize them: positions and forces are
reached through `__restrict__`-qualified typed pointers, and force accumulation
moves out of the innermost loop into scalar locals. Redundant `pow()` calls are
removed from the `nm/cut`, `lj/pirani`, and `mie/cut` families.

## Contents

```
src/                       modified and new pair styles, laid out by LAMMPS package
  OPENMP/                    36 new + 4 updated pair_*_omp.{cpp,h}
  <PACKAGE>/                 225 rewritten pair_*.cpp across 38 packages
scripts/apply_patch.sh     applies the package to an existing LAMMPS tree
tools/verify/              decks and harness for the 11 styles the upstream
                           regression suite does not reach
tools/bench/               benchmark inputs and the measurement harness
Paper/tableS6_a3_files.csv file-level inventory: path, package, pair_style name,
                           and which code markers each file carries
INSTALL.md                 build instructions
LICENSE                    GPLv2
```

## Installing

```bash
git clone -b develop https://github.com/lammps/lammps.git
cd lammps && git checkout 91d4111a9        # the commit this was developed against
cd ..
bash lammps-CO/scripts/apply_patch.sh ./lammps

cd lammps && mkdir build && cd build
cmake ../cmake -D CMAKE_BUILD_TYPE=Release -D PKG_OPENMP=on \
      -D BUILD_OMP=on -D BUILD_MPI=on
cmake --build . -j
```

`apply_patch.sh` walks this repository's `src/` and replaces matching files in
the target tree. Packages absent from the target are reported, not created. No
CMake changes are needed — LAMMPS discovers `_omp.cpp` files automatically.

## Verifying

Build the same LAMMPS commit twice, once unmodified and once with this package
applied, then run both through LAMMPS's own regression suite:

```bash
cmake ../cmake -D ENABLE_TESTING=on ...     # in each build directory
ctest -R "PairStyle"
```

The two builds must return identical outcomes on all 325 pair-style tests.
For the eleven styles that suite does not reach:

```bash
bash tools/verify/compare_decks.sh tools/verify/decks.tsv
```

which runs each deck under both binaries and compares thermo output character
for character. See `tools/verify/README.md`.

## Known limitations

Three files were rewritten and then reverted, because the rewrite did not
preserve floating-point operation order for them: `pair_buck_coul_cut`,
`pair_kolmogorov_crespi_z`, and `pair_polymorphic`. They ship unmodified.

`pair_kolmogorov_crespi_z.cpp` carries one unrelated one-line change: stock
LAMMPS aborts on `pair_style hybrid/overlay kolmogorov/crespi/z` under
`-sf omp`, because `hybrid/overlay` is registered under two names and the
sub-style compares `force->pair_style` for exact equality. The fix replaces that
with a prefix match and is offered upstream separately.
