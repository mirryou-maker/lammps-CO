# LLM-Assisted LAMMPS OpenMP Optimization

**Repository for**: "AI-Assisted Systematic Optimization of the LAMMPS Molecular Dynamics Simulator via Large Language Model Coding Agents"  
**Author**: Chun-Yeol You, DGIST, Republic of Korea  
**Contact**: cyyou@dgist.ac.kr

---

## Contents

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
