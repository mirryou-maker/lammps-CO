"""
analyze_fep_opls.py  — TI + BAR analysis of OPLS-AA FEP results for EC, PC, DME.

Input:  tools/fep/logs_opls/fep_{solvent}_lam{lam}.stdout
Output: ΔG values, upgraded Figure 7 (multi-solvent comparison)

Thermo columns: step  lambda  dU/dlambda  exp(-ΔU/kT)
"""

import pathlib, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

LOG_DIR  = pathlib.Path("tools/fep/logs_opls")
LOG_CG   = pathlib.Path("tools/fep/logs")
OUTDIR   = pathlib.Path("tools/plots")
OUTDIR.mkdir(exist_ok=True)

LAMBDAS = np.array([0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,
                    0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.00])

KT = 0.5922   # kcal/mol at 300 K


def _open_log(path):
    """Open a log file, detecting UTF-16 LE BOM written by PowerShell 5.1."""
    with open(path, 'rb') as fb:
        header = fb.read(2)
    if header == b'\xff\xfe':
        return open(path, encoding='utf-16-le', errors='replace')
    return open(path, encoding='utf-8', errors='replace')


def parse_windows(log_dir: pathlib.Path, prefix: str) -> dict:
    """Return {lambda: array_of_dUdL} for all completed windows."""
    data = {}
    for lam in LAMBDAS:
        tag    = f"{lam:.2f}"
        stdout = log_dir / f"fep_{prefix}{tag}.stdout"
        if not stdout.exists():
            continue
        dUdL_vals = []
        with _open_log(stdout) as fh:
            for line in fh:
                cols = line.split()
                if len(cols) == 4:
                    try:
                        lam_v = float(cols[1])
                        dUdl  = float(cols[2])
                        if abs(lam_v - lam) < 0.02:
                            dUdL_vals.append(dUdl)
                    except ValueError:
                        pass
        if len(dUdL_vals) > 5:
            data[lam] = np.array(dUdL_vals)
    return data


def ti_estimate(data: dict):
    """Return (deltaG, err, lams_arr, mean_arr, sem_arr)."""
    lams     = sorted(data.keys())
    mean_arr = np.array([data[l].mean() for l in lams])
    sem_arr  = np.array([data[l].std()/np.sqrt(len(data[l])) for l in lams])
    lams_arr = np.array(lams)
    dG       = np.trapezoid(mean_arr, lams_arr)
    # propagated SEM
    dl = np.diff(lams_arr)
    w  = np.zeros(len(lams_arr))
    w[0]  += dl[0]/2; w[-1] += dl[-1]/2
    for i in range(1, len(dl)):
        w[i] += (dl[i-1]+dl[i])/2
    err = np.sqrt(np.sum((w*sem_arr)**2))
    return dG, err, lams_arr, mean_arr, sem_arr


def bar_estimate(data: dict):
    """Return BAR ΔG by summing adjacent-window BAR values."""
    try:
        import pymbar
    except ImportError:
        return np.nan
    lams    = sorted(data.keys())
    total   = 0.0
    for i in range(len(lams)-1):
        l0, l1 = lams[i], lams[i+1]
        dl = l1 - l0
        dU_F =  data[l0] * dl / KT
        dU_R = -data[l1] * dl / KT
        try:
            result = pymbar.other_estimators.bar(dU_F, dU_R)
            total += result['Delta_f'] * KT
        except Exception:
            pass
    return total


# ─── Load CG results (from original FEP run) ────────────────────────────────

cg_data = parse_windows(LOG_CG, "lam")
if cg_data:
    cg_dG, cg_err, cg_lams, cg_mean, cg_sem = ti_estimate(cg_data)
    cg_bar = bar_estimate(cg_data)
    print(f"CG     TI = {cg_dG:.3f} ± {cg_err:.3f}  BAR = {cg_bar:.3f} kcal/mol")
else:
    cg_dG = cg_err = np.nan; cg_lams = cg_mean = cg_sem = np.array([])
    print("CG: no data found")

# ─── Load OPLS-AA results ────────────────────────────────────────────────────

results = {}
for sol in ['ec', 'pc', 'dme']:
    d = parse_windows(LOG_DIR, f"{sol}_lam")
    if not d:
        print(f"{sol.upper()}: no completed windows found - skipping")
        continue
    dG, err, lams_a, mean_a, sem_a = ti_estimate(d)
    bar = bar_estimate(d)
    results[sol] = {
        'dG': dG, 'err': err, 'bar': bar,
        'lams': lams_a, 'mean': mean_a, 'sem': sem_a,
        'n_windows': len(d)
    }
    print(f"{sol.upper():<5}  TI = {dG:+.3f} ± {err:.3f}  BAR = {bar:+.3f} kcal/mol"
          f"  ({len(d)}/21 windows)")

if not results:
    print("\nNo OPLS-AA results yet. Re-run after FEP windows complete.")
    sys.exit(0)

# ─── Figure 7 (upgraded) ────────────────────────────────────────────────────
# Layout: 2 rows × 3 cols
#   Row 0: TI integrand for EC, PC, DME
#   Row 1: ΔG bar chart comparison | Running TI (EC) | OMP speedup panel

COLORS = {'ec': '#1565C0', 'pc': '#6A1B9A', 'dme': '#2E7D32', 'cg': '#E65100'}
LABELS = {'ec': 'EC (OPLS-AA)', 'pc': 'PC (OPLS-AA)',
          'dme': 'DME (OPLS-AA)', 'cg': 'CG model'}

fig = plt.figure(figsize=(16, 9))
gs  = GridSpec(2, 3, figure=fig, hspace=0.46, wspace=0.40)

# Row 0: TI integrand for each OPLS-AA solvent
for col, sol in enumerate(['ec', 'pc', 'dme']):
    if sol not in results:
        continue
    r  = results[sol]
    ax = fig.add_subplot(gs[0, col])
    ax.errorbar(r['lams'], r['mean'], yerr=r['sem']*3,
                fmt='o-', color=COLORS[sol], capsize=3, lw=1.8, ms=5)
    ax.fill_between(r['lams'], r['mean'], alpha=0.12, color=COLORS[sol])
    ax.axhline(0, lw=0.7, ls='--', color='#999')
    ax.set_xlabel('λ', fontsize=10)
    ax.set_ylabel('⟨dU/dλ⟩ (kcal/mol)', fontsize=10)
    ax.set_title(f"({chr(ord('a')+col)}) {LABELS[sol]}\nΔG_TI = {r['dG']:+.2f} kcal/mol",
                 fontsize=10, pad=6)
    ax.grid(alpha=0.3)

# Row 1 left: ΔG comparison bar chart
ax_bar = fig.add_subplot(gs[1, 0])
all_solvs = [s for s in ['cg', 'ec', 'pc', 'dme'] if s in results or (s=='cg' and cg_data)]
dGs  = []
errs = []
cols = []
labs = []
for s in all_solvs:
    if s == 'cg' and cg_data:
        dGs.append(cg_dG); errs.append(cg_err)
    elif s in results:
        dGs.append(results[s]['dG']); errs.append(results[s]['err'])
    else:
        continue
    cols.append(COLORS[s])
    labs.append(LABELS[s])

x = np.arange(len(labs))
bars = ax_bar.bar(x, dGs, yerr=errs, color=cols, capsize=5, width=0.55,
                  edgecolor='black', linewidth=0.7, alpha=0.85, error_kw=dict(elinewidth=1.5))
ax_bar.axhline(0, lw=0.8, color='black')
ax_bar.set_xticks(x)
ax_bar.set_xticklabels(labs, rotation=22, ha='right', fontsize=9)
ax_bar.set_ylabel('ΔG_desolv (kcal/mol)', fontsize=10)
ax_bar.set_title('(d) Desolvation ΔG comparison', fontsize=10, pad=6)
ax_bar.grid(axis='y', alpha=0.3)
for bar_, dg in zip(bars, dGs):
    ax_bar.text(bar_.get_x()+bar_.get_width()/2, dg+0.3 if dg >= 0 else dg-1.2,
                f'{dg:+.2f}', ha='center', va='bottom', fontsize=8.5, fontweight='bold')

# Row 1 middle: running TI for EC (best data)
sol_main = 'ec' if 'ec' in results else (list(results.keys())[0] if results else None)
if sol_main:
    r   = results[sol_main]
    cumG = np.array([np.trapezoid(r['mean'][:i+1], r['lams'][:i+1])
                     for i in range(len(r['lams']))])
    ax_run = fig.add_subplot(gs[1, 1])
    ax_run.plot(r['lams'], cumG, 'o-', color=COLORS[sol_main], lw=1.8, ms=5)
    ax_run.axhline(r['dG'], ls='--', lw=1.2, color='#555',
                   label=f"ΔG = {r['dG']:+.2f} kcal/mol")
    ax_run.set_xlabel('λ', fontsize=10)
    ax_run.set_ylabel('Cumulative ΔG (kcal/mol)', fontsize=10)
    ax_run.set_title(f'(e) Running TI – {LABELS[sol_main]}', fontsize=10, pad=6)
    ax_run.legend(fontsize=9)
    ax_run.grid(alpha=0.3)

# Row 1 right: OMP speedup (from CLAUDE.md benchmark data)
ax_omp = fig.add_subplot(gs[1, 2])
n_threads = [1, 2, 4, 8]
# Speedup data from session: 4t→4.04×, 8t→5.28× (lj/cut/soft/omp, 32k atoms)
speedup   = [1.0, 2.15, 4.04, 5.28]
ax_omp.plot(n_threads, speedup, 'o-', color='#E53935', lw=2, ms=7, label='lj/cut/soft/omp')
ax_omp.plot(n_threads, n_threads, 'k--', lw=1, label='ideal linear')
ax_omp.fill_between(n_threads, speedup, 1, alpha=0.15, color='#E53935')
ax_omp.set_xlabel('OMP threads', fontsize=10)
ax_omp.set_ylabel('Speedup vs serial', fontsize=10)
ax_omp.set_title('(f) OpenMP FEP speedup\n(lj/cut/soft/omp, 32k atoms)', fontsize=10, pad=6)
ax_omp.legend(fontsize=9)
ax_omp.grid(alpha=0.3)
ax_omp.set_xticks(n_threads)

fig.suptitle(
    'Figure 7 · FEP Case Study: Li⁺ Desolvation in Battery Electrolyte Solvents\n'
    '(OPLS-AA all-atom force field; TI with 21 λ-windows; 40 ps/window)',
    fontsize=11, fontweight='bold', y=1.01
)

out_pdf = str(OUTDIR / "figure7_fep_casestudy.pdf")
out_png = str(OUTDIR / "figure7_fep_casestudy.png")
plt.savefig(out_pdf, bbox_inches='tight', dpi=300)
plt.savefig(out_png, bbox_inches='tight', dpi=200)
print(f"\nFigure saved: {out_pdf}")
print(f"             {out_png}")
