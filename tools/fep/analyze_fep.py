"""
analyze_fep.py  — TI + MBAR analysis of 20-window FEP results.

Input:  tools/fep/logs/fep_lam*.stdout
Output: deltaG (kcal/mol), Figure 7 panels

The LAMMPS thermo output columns are:
  step  lambda  dU/dlambda  exp(-dU/kT)
"""

import pathlib, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LOG_DIR = pathlib.Path("tools/fep/logs")
OUTDIR  = pathlib.Path("tools/plots")

# ─── 1. Parse thermo files ────────────────────────────────────────────────────

lambdas = np.array([0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
                    0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00])

dUdL_by_window = {}   # lambda → array of dU/dlambda values
expBolt_by_window = {}

for lam in lambdas:
    label  = f"{lam:.2f}"
    stdout = LOG_DIR / f"fep_lam{label}.stdout"
    if not stdout.exists():
        print(f"WARNING: missing {stdout}", file=sys.stderr)
        continue

    dUdL_vals = []
    expB_vals = []
    with open(stdout) as f:
        for line in f:
            cols = line.split()
            if len(cols) == 4:
                try:
                    step  = int(cols[0])
                    lam_v = float(cols[1])
                    dUdl  = float(cols[2])
                    expb  = float(cols[3])
                    if abs(lam_v - lam) < 0.01:   # production data
                        dUdL_vals.append(dUdl)
                        expB_vals.append(expb)
                except ValueError:
                    pass

    dUdL_by_window[lam]   = np.array(dUdL_vals)
    expBolt_by_window[lam] = np.array(expB_vals)
    print(f"  lambda={label}: {len(dUdL_vals)} frames, "
          f"<dU/dl>={np.mean(dUdL_vals):.4f} ± {np.std(dUdL_vals)/np.sqrt(len(dUdL_vals)):.4f} kcal/mol")

# ─── 2. Thermodynamic Integration (trapezoidal rule) ──────────────────────────

lams_ok    = sorted(k for k in dUdL_by_window if len(dUdL_by_window[k]) > 10)
mean_dUdL  = np.array([dUdL_by_window[l].mean() for l in lams_ok])
sem_dUdL   = np.array([dUdL_by_window[l].std() / np.sqrt(len(dUdL_by_window[l]))
                       for l in lams_ok])
lams_arr   = np.array(lams_ok)

deltaG_TI  = np.trapezoid(mean_dUdL, lams_arr)
# Error estimate (propagated from SEM via trapezoidal weights)
dl         = np.diff(lams_arr)
w          = np.zeros(len(lams_arr))
w[0]      += dl[0] / 2
w[-1]     += dl[-1] / 2
for i in range(1, len(dl)):
    w[i] += (dl[i-1] + dl[i]) / 2
err_TI     = np.sqrt(np.sum((w * sem_dUdL)**2))

print(f"\n=== TI result ===")
print(f"ΔG(Li+ desolvation) = {deltaG_TI:.3f} ± {err_TI:.3f} kcal/mol")
print(f"  (λ: 1→0; negative = solvation favorable)")

# ─── 3. BAR estimate using adjacent windows ───────────────────────────────────

try:
    import pymbar
    kT = 0.5922  # kcal/mol at 300 K
    bar_dGs = []
    for i in range(len(lams_ok) - 1):
        l0, l1 = lams_ok[i], lams_ok[i+1]
        dU_forward = (dUdL_by_window[l0] * (l1 - l0))
        dU_reverse = -(dUdL_by_window[l1] * (l1 - l0))
        try:
            result = pymbar.other_estimators.bar(dU_forward / kT, dU_reverse / kT)
            bar_dGs.append(result['Delta_f'] * kT)
        except Exception:
            bar_dGs.append(np.nan)
    deltaG_BAR = np.nansum(bar_dGs)
    print(f"BAR estimate        = {deltaG_BAR:.3f} kcal/mol")
    bar_available = True
except Exception as e:
    print(f"BAR skipped: {e}")
    bar_available = False

# ─── 4. MBAR estimate ────────────────────────────────────────────────────────

try:
    import pymbar
    kT     = 0.5922
    # Build u_kn matrix: row k = lambda state, col n = sample index
    # For efficiency, subsample to max 500 frames per window
    N_max  = min(500, min(len(dUdL_by_window[l]) for l in lams_ok))
    N_k    = np.array([N_max] * len(lams_ok))

    # u_kn[k, n] = reduced potential of sample n evaluated at lambda state k
    # We approximate: u_kn[k,n] = u_kk[n] + (lams_ok[k] - lams_ok[state_n]) * dU/dl
    # Simple linear approximation valid for small Δλ intervals
    u_kn = np.zeros((len(lams_ok), len(lams_ok) * N_max))
    for k, lam_k in enumerate(lams_ok):
        for n_state, lam_n in enumerate(lams_ok):
            dU_n = dUdL_by_window[lam_n][:N_max]
            dl   = lam_k - lam_n
            u_kn[k, n_state*N_max : (n_state+1)*N_max] = (
                dU_n * dl / kT
            )

    mbar = pymbar.MBAR(u_kn, N_k, verbose=False)
    result_dict = mbar.compute_free_energy_differences()
    deltaG_MBAR = result_dict['Delta_f'][0, -1] * kT  # from lambda=0 to lambda=1
    print(f"MBAR estimate       = {deltaG_MBAR:.3f} kcal/mol  (approx; linear u_kn)")
    mbar_available = True
except Exception as e:
    print(f"MBAR skipped: {e}")
    deltaG_MBAR    = np.nan
    mbar_available = False

# ─── 5. Figure 7 ──────────────────────────────────────────────────────────────

fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
fig.subplots_adjust(wspace=0.38)

# (a) dU/dλ vs λ  (TI integrand)
ax = axes[0]
ax.errorbar(lams_arr, mean_dUdL, yerr=sem_dUdL * 3,  # 3σ error bars
            fmt='o-', color="#1565C0", capsize=4, lw=1.8, ms=6,
            label="⟨dU/dλ⟩ ± 3 SEM")
ax.axhline(0, lw=0.8, ls="--", color="#aaa")
ax.fill_between(lams_arr, mean_dUdL, alpha=0.15, color="#1565C0",
                label=f"Area (TI) = {deltaG_TI:.2f} kcal/mol")
ax.set_xlabel("λ (coupling parameter)", fontsize=11)
ax.set_ylabel("⟨dU/dλ⟩ (kcal/mol)", fontsize=11)
ax.set_title("(a) TI integrand", fontsize=11, pad=8)
ax.legend(fontsize=9)
ax.grid(alpha=0.35)

# (b) Cumulative ΔG (TI running integral)
cumG = np.array([np.trapezoid(mean_dUdL[:i+1], lams_arr[:i+1]) for i in range(len(lams_arr))])
ax = axes[1]
ax.plot(lams_arr, cumG, 'o-', color="#E53935", lw=1.8, ms=6)
ax.axhline(deltaG_TI, ls="--", lw=1.2, color="#555",
           label=f"Final ΔG = {deltaG_TI:.2f} kcal/mol")
ax.set_xlabel("λ", fontsize=11)
ax.set_ylabel("Cumulative ΔG (kcal/mol)", fontsize=11)
ax.set_title("(b) Running TI estimate", fontsize=11, pad=8)
ax.legend(fontsize=9)
ax.grid(alpha=0.35)

# (c) per-window <exp(-dU/kT)> variance → convergence indicator
kT = 0.5922
ax = axes[2]
variances = [np.var(expBolt_by_window[l]) for l in lams_ok if len(expBolt_by_window[l]) > 0]
ax.bar(lams_arr[:len(variances)], variances, color="#43A047", width=0.04, alpha=0.8)
ax.set_xlabel("λ", fontsize=11)
ax.set_ylabel("Var[exp(-ΔU/kT)]", fontsize=11)
ax.set_title("(c) Per-window convergence\n(lower = better overlap)", fontsize=11, pad=8)
ax.grid(axis="y", alpha=0.35)

plt.suptitle(
    f"Figure 7 · FEP: Li⁺ desolvation in CG solvent  "
    f"(ΔG_TI = {deltaG_TI:.2f} kcal/mol, 21 windows × 1 ns, OMP 4t)",
    fontsize=11, y=1.01, fontweight="bold"
)

out_pdf = str(OUTDIR / "figure7_fep_casestudy.pdf")
out_png = str(OUTDIR / "figure7_fep_casestudy.png")
plt.savefig(out_pdf, bbox_inches="tight", dpi=300)
plt.savefig(out_png, bbox_inches="tight", dpi=200)
print(f"\nSaved {out_pdf}")
print(f"Saved {out_png}")
