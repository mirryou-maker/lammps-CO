# LLM-Assisted LAMMPS OpenMP Optimization

**Repository for**: "AI-Assisted Systematic Optimization of the LAMMPS Molecular Dynamics Simulator via Large Language Model Coding Agents"  
**Author**: Chun-Yeol You, DGIST, Republic of Korea  
**Contact**: cyyou@dgist.ac.kr

---

## Contents

### `src/OPENMP/` — New OpenMP pair_style variants (A-1, this work)
36 new `pair_*_omp.{h,cpp}` file pairs created by LLM-assisted code generation.
Copy into `lammps/src/OPENMP/` and rebuild with `-DPKG_OPENMP=on` to enable.

| Style group | Files |
|---|---|
| born/coul (dsf, dsf_cs, long_cs, wolf_cs) | 4 pairs |
| bpm_spring, buck_coul_long_cs | 2 pairs |
| coul (ctip, cut_dielectric, cut_soft_gapsys, esp, exclude, long_cs, long_dielectric, slater_cut, slater_long, tt, wolf_cs) | 11 pairs |
| dpd_fdt, kolmogorov_crespi_z, lebedeva_z | 3 pairs |
| lj_charmmfsw (coul_charmmfsh, coul_long) | 2 pairs |
| lj_class2 (coul_cut_soft, coul_long_cs, coul_long_soft, soft) | 4 pairs |
| lj_cut (coul_esp, coul_long_cs, coul_msm_dielectric, dipole_long) | 4 pairs |
| lj_expand_coul_long, mdpd, morse_soft, nm_cut_split, rheo_solid, thole | 6 pairs |

### `src/` — A-3 optimized standard pair files (15 files)
Standard `pair_*.cpp` files backported with `__restrict__`/`_noalias` hints and
`fxtmp`/`fytmp`/`fztmp` register-accumulator pattern to enable compiler auto-vectorization.
Copy into `lammps/src/` (replaces originals — bit-identical results verified).

Files: `pair_born`, `pair_buck`, `pair_buck_coul_cut`, `pair_coul_cut`, `pair_coul_debye`,
`pair_coul_dsf`, `pair_coul_wolf`, `pair_lj_cut`, `pair_lj_cut_coul_cut`,
`pair_lj_expand`, `pair_morse`, `pair_soft`, `pair_table`, `pair_yukawa`, `pair_zbl`

---

### `tools/bench/` — Performance benchmark harness (PowerShell)
- `Run-Benchmark.ps1` — Runs LAMMPS N times, parses Loop/Pair/Neigh/Comm times, saves JSON + CSV
- `LammpsLog.psm1` — LAMMPS stdout parser module
- `Build-Lammps.ps1` — CMake + MSVC build helper
- `results/` — Benchmark JSON results and `summary.csv`

### `tools/accel/` — Accelerator coverage checker
- `Check-AccelCoverage.ps1` — Scans a LAMMPS input script and reports which pair/bond/angle styles have unused OMP/GPU/KOKKOS variants

### `tools/fep/` — Free Energy Perturbation workflow
- `run_fep_opls.ps1` — 40 ps/window FEP runner (EC, PC, DME — 63 windows total)
- `run_fep_opls_long.ps1` — 200 ps/window extended runner
- `analyze_fep_opls.py` — TI + BAR estimator using pymbar
- `data/` — LAMMPS data files for CG and OPLS-AA Li⁺/solvent systems
  - `lj_ec_li.data` — CG model (300 EC beads + 1 Li⁺, 301 atoms)
  - `ec_li_aa.data`, `pc_li_aa.data`, `dme_li_aa.data` — OPLS-AA all-atom systems

### `tools/plots/` — Figure generation scripts (matplotlib)
- `fig1_workflow.py` — LLM-assisted optimization pipeline
- `fig2_omp_coverage.py` — OMP pair_style coverage (51% → 63%)
- `fig3_build_flags.py` — Compiler flag benchmark results
- `fig4_omp_scaling.py` — OpenMP scaling and Amdahl's law fit
- `fig5_a3_backport.py` — A-3 restrict/_noalias optimization
- `fig6_summary.py` — Combined summary (waterfall, tornado, heatmap, donut)
- `fig7_fep_casestudy.py` — FEP case study (TI integrands, ΔG bar, convergence)
- `run_all_figures.py` — Runs fig1–fig6 in sequence

---

## Key Results

| Optimization | Speedup |
|---|---|
| OMP pair coverage | 157 → 192/307 styles (+35, +11.4 pp) |
| OMP 4 threads (32k atoms) | 4.31× |
| OMP 8 threads (32k atoms) | 4.93× |
| AVX2 + /fp:fast build flag | +6.3% |
| A-3 restrict backport (lj/cut) | +4.2% |
| FEP OMP 4 threads (lj/cut/soft) | 4.04× |
| FEP OMP 8 threads (lj/cut/soft) | 5.28× |

### Li⁺ Desolvation Free Energies (OPLS-AA, 200 ps/window)

| Solvent | ΔG_TI (kcal/mol) | ΔG_BAR (kcal/mol) |
|---|---|---|
| Ethylene carbonate (EC) | −183.7 ± 0.1 | −168.4 |
| Propylene carbonate (PC) | −182.1 ± 0.1 | −167.2 |
| 1,2-Dimethoxyethane (DME) | −159.4 ± 0.1 | −147.0 |

---

## Requirements

- LAMMPS (develop branch, OPENMP package enabled)
- Python 3.x with: `numpy`, `matplotlib`, `scipy`, `pymbar`
- PowerShell 5.1+ (Windows) or pwsh (Linux/Mac)
- CMake 3.16+, MSVC 2022 (Windows) or GCC/Clang (Linux/Mac)
