# Cover Letter

**Chun-Yeol You**
Department of Physics and Chemistry
DGIST (Daegu Gyeongbuk Institute of Science and Technology)
Daegu, Republic of Korea
cyyou@dgist.ac.kr

---

**To the Editors of *Computer Physics Communications***

Dear Editors,

I submit for your consideration the manuscript **"LAMMPS-CO: extended OpenMP coverage and compiler-friendly force kernels for 261 LAMMPS pair styles"** as a Computer Programs in Physics (CPiP) paper.

## What the program does

LAMMPS parallelizes across MPI ranks by spatial decomposition, but within a rank a pair style is threaded only if an OpenMP variant of that style exists — and fewer than half of the 307 pair styles have one. On a many-core node, a simulation using one of the remaining styles leaves most of the node idle. Separately, the serial force kernels are written in a form the compiler cannot vectorize: positions and forces are reached through unqualified `double **` pointers, so aliasing cannot be ruled out, and the per-atom force array is updated inside the innermost loop, creating a loop-carried dependence.

LAMMPS-CO addresses both. It supplies 36 new OpenMP pair-style variants, raising coverage from 157 to 193 of 307 styles, and rewrites the force kernels of 225 further pair files across 38 packages to use `__restrict__`-qualified typed pointers and local accumulators. Redundant `pow()` evaluations are removed from the `nm/cut`, `lj/pirani`, and `mie/cut` families. Measured gains are 4.31× at four OpenMP threads and 4.93× at eight on a 32,000-atom benchmark, 4–8% on pair time from the kernel rewrite, and 19.9–40.8% within the styles where `pow()` was eliminated.

## Why it belongs in CPC

The program is a direct, reusable contribution to the most widely used open-source molecular dynamics engine, and it targets a gap that affects every LAMMPS user running a style without an accelerated variant. It also illustrates a programming technique of current interest to the computational physics community: the code was generated with the assistance of a large language model coding agent, working from the reference implementations already present in the OPENMP and OPT packages. We report that protocol, its measured error rate, and — importantly — the verification that caught its errors.

That verification is the part I would draw your attention to. Rather than assert numerical equivalence, we build the unmodified and modified trees from the same LAMMPS commit and run both through LAMMPS's own `unittest/force-styles` regression suite: 325 pair-style tests, identical outcomes. Eleven modified styles are not reachable by that suite; we distribute dedicated decks and a comparison harness for those, so every modified style is checked. The process found ten files that did not compile in packages absent from our original build, one file whose force expression had a multiplication transcribed as a division, and three files where the rewrite did not preserve operation order — the last of which we reverted rather than ship. All of this is reported in the manuscript. We take the view that a software paper is more useful for saying what the verification found than for claiming it found nothing.

The exercise also surfaced a defect in unmodified LAMMPS: `pair_style hybrid/overlay kolmogorov/crespi/z` aborts under `-sf omp`, because `hybrid/overlay` is registered under two names and the sub-style tests `force->pair_style` for exact equality. The upstream test suite cannot see this, as the OpenMP branch is skipped when no OMP variant of the sub-style exists. We have submitted the fix upstream as LAMMPS pull request #5099.

## Availability

The code is public at <https://github.com/mirryou-maker/lammps-CO> under GPLv2, matching LAMMPS, and is being deposited in the CPC Program Library on Mendeley Data. It applies to an existing LAMMPS tree with a single script and needs no build-system changes. The complete file-level inventory of modified styles, the verification decks, and the benchmark harness are distributed with it.

## Statements

This manuscript is original, is not under consideration elsewhere, and has not been published previously. An earlier version was submitted to *Archives of Computational Methods in Engineering* and declined without review as outside that journal's scope (it publishes state-of-the-art reviews); it has not been peer reviewed. The author declares no competing interests.

Thank you for considering this submission.

Sincerely,

Chun-Yeol You
DGIST
