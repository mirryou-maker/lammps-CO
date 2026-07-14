"""
Figure 7: FEP Case Study — Li+ desolvation in battery electrolyte solvents.
Uses 200 ps/window extended run data (logs_opls_long/).
(a-c) TI integrands for EC, PC, DME
(d)   ΔG bar chart (TI vs BAR, 40 ps vs 200 ps)
(e)   Running TI for EC
(f)   OMP speedup for lj/cut/soft/omp
"""

import sys, pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

plt.rcParams.update({
    'font.size': 13,
    'axes.titlesize': 14,
    'axes.labelsize': 13,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.titlesize': 15,
    'lines.linewidth': 2.2,
    'lines.markersize': 7,
})

# ── paths ──────────────────────────────────────────────────────────────────
ROOT    = pathlib.Path(__file__).parent.parent.parent   # D:/Claude-Code-R/LAMMPS-CO
LOG_LONG = ROOT / "tools/fep/logs_opls_long"
LOG_SHORT = ROOT / "tools/fep/logs_opls"
OUTDIR   = ROOT / "tools/plots"
OUTDIR.mkdir(exist_ok=True)

LAMBDAS = np.array([0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,
                    0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.00])
KT = 0.5922

COLORS = {'ec': '#1565C0', 'pc': '#6A1B9A', 'dme': '#2E7D32'}
SOLVENT_LABELS = {'ec': 'EC (ε=89.8)', 'pc': 'PC (ε=64.9)', 'dme': 'DME (ε=7.2)'}

# ── helpers ────────────────────────────────────────────────────────────────
def _open_log(path):
    with open(path, 'rb') as fb:
        header = fb.read(2)
    if header == b'\xff\xfe':
        return open(path, encoding='utf-16-le', errors='replace')
    return open(path, encoding='utf-8', errors='replace')


def parse_windows(log_dir, prefix):
    data = {}
    for lam in LAMBDAS:
        fp = log_dir / f"fep_{prefix}{lam:.2f}.stdout"
        if not fp.exists():
            continue
        vals = []
        with _open_log(fp) as fh:
            for line in fh:
                cols = line.split()
                if len(cols) == 4:
                    try:
                        if abs(float(cols[1]) - lam) < 0.02:
                            vals.append(float(cols[2]))
                    except ValueError:
                        pass
        if len(vals) > 5:
            data[lam] = np.array(vals)
    return data


def ti_bar(data):
    lams     = np.array(sorted(data.keys()))
    means    = np.array([data[l].mean() for l in lams])
    sems     = np.array([data[l].std() / np.sqrt(len(data[l])) for l in lams])
    dG       = np.trapezoid(means, lams)
    dl       = np.diff(lams)
    w        = np.zeros(len(lams))
    w[0]    += dl[0] / 2;  w[-1] += dl[-1] / 2
    for i in range(1, len(dl)):
        w[i] += (dl[i-1] + dl[i]) / 2
    err = np.sqrt(np.sum((w * sems) ** 2))
    try:
        import pymbar
        total = 0.0
        for i in range(len(lams) - 1):
            dl2 = lams[i+1] - lams[i]
            try:
                total += pymbar.other_estimators.bar(
                    data[lams[i]]   * dl2 / KT,
                    -data[lams[i+1]] * dl2 / KT)['Delta_f'] * KT
            except Exception:
                pass
        bar = total
    except ImportError:
        bar = float('nan')
    return dG, err, lams, means, sems, bar


# ── load 200 ps results ─────────────────────────────────────────────────────
results = {}
for sol in ['ec', 'pc', 'dme']:
    d = parse_windows(LOG_LONG, f"{sol}_lam")
    if not d:
        print(f"  {sol.upper()}: no data in {LOG_LONG} — trying short run")
        d = parse_windows(LOG_SHORT, f"{sol}_lam")
    if d:
        dG, err, lams, means, sems, bar = ti_bar(d)
        results[sol] = dict(dG=dG, err=err, lams=lams, means=means, sems=sems,
                            bar=bar, n=len(d))
        print(f"  {sol.upper()}: TI={dG:+.2f}±{err:.2f}  BAR={bar:+.2f}  ({len(d)}/21 win)")

if not results:
    print("No data found. Exiting.")
    sys.exit(1)

# ── confirmed 40 ps values (from short run, for comparison panel) ──────────
SHORT_TI  = {'ec': -184.703, 'pc': -182.848, 'dme': -159.877}
SHORT_BAR = {'ec': -169.116, 'pc': -167.820, 'dme': -147.434}
SHORT_ERR = {'ec':  0.276,   'pc':  0.286,   'dme':  0.275  }

# ── figure ─────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(22, 13))
gs  = GridSpec(2, 3, figure=fig, hspace=0.58, wspace=0.48)

# Row 0: TI integrands ──────────────────────────────────────────────────────
for col, sol in enumerate(['ec', 'pc', 'dme']):
    if sol not in results:
        continue
    r  = results[sol]
    ax = fig.add_subplot(gs[0, col])
    ax.errorbar(r['lams'], r['means'], yerr=r['sems'] * 3,
                fmt='o-', color=COLORS[sol], capsize=4, lw=2.2, ms=7,
                label='⟨dU/dλ⟩ ±3σ')
    ax.fill_between(r['lams'], r['means'], alpha=0.13, color=COLORS[sol])
    ax.axhline(0, lw=1.0, ls='--', color='#999')

    # annotate the prominent peak
    peak_idx = np.argmin(r['means'])
    ax.annotate(f"min: {r['means'][peak_idx]:.0f} kcal/mol",
                xy=(r['lams'][peak_idx], r['means'][peak_idx]),
                xytext=(r['lams'][peak_idx] + 0.12, r['means'][peak_idx] + 40),
                fontsize=10.5, color=COLORS[sol],
                arrowprops=dict(arrowstyle='->', color=COLORS[sol], lw=1.5))

    ax.set_xlabel('λ', fontsize=13)
    ax.set_ylabel('⟨dU/dλ⟩ (kcal/mol)', fontsize=13)
    ax.set_title(
        f"({'abc'[col]}) {SOLVENT_LABELS[sol]}\n"
        f"ΔG_TI = {r['dG']:+.1f} kcal/mol  ({r['n']}/21 win)",
        fontsize=14, pad=8)
    ax.legend(fontsize=11, loc='lower right')
    ax.grid(alpha=0.3)

# Row 1 left: ΔG comparison (200 ps TI vs BAR, 40 ps TI for comparison) ───
ax_bar = fig.add_subplot(gs[1, 0])
solvents_list = [s for s in ['ec', 'pc', 'dme'] if s in results]
x = np.arange(len(solvents_list))
bar_w = 0.26

b1 = ax_bar.bar(x - bar_w, [SHORT_TI[s]  for s in solvents_list],
                bar_w, color='#90A4AE', label='TI 40 ps', zorder=3,
                yerr=[SHORT_ERR[s] for s in solvents_list], capsize=4,
                error_kw=dict(elinewidth=1.5))
b2 = ax_bar.bar(x,          [results[s]['dG']  for s in solvents_list],
                bar_w, color='#1565C0', label='TI 200 ps', zorder=3,
                yerr=[results[s]['err'] for s in solvents_list], capsize=4,
                error_kw=dict(elinewidth=1.5))
b3 = ax_bar.bar(x + bar_w,  [results[s]['bar'] for s in solvents_list],
                bar_w, color='#E53935', label='BAR 200 ps', zorder=3)

ax_bar.axhline(0, lw=0.8, color='black')
ax_bar.set_xticks(x)
ax_bar.set_xticklabels([SOLVENT_LABELS[s] for s in solvents_list],
                       rotation=18, ha='right', fontsize=12)
ax_bar.set_ylabel('ΔG_desolv (kcal/mol)', fontsize=13)
ax_bar.set_title('(d) Li⁺ desolvation free energy\nTI vs BAR, 40 ps vs 200 ps',
                 fontsize=14, pad=8)
ax_bar.legend(fontsize=11, loc='lower right')
ax_bar.grid(axis='y', alpha=0.3)

# Row 1 middle: running TI for EC ──────────────────────────────────────────
sol_main = 'ec' if 'ec' in results else list(results.keys())[0]
r = results[sol_main]
cumG = np.array([np.trapezoid(r['means'][:i+1], r['lams'][:i+1])
                 for i in range(len(r['lams']))])
ax_run = fig.add_subplot(gs[1, 1])
ax_run.plot(r['lams'], cumG, 'o-', color=COLORS[sol_main], lw=2.2, ms=7)
ax_run.axhline(r['dG'], ls='--', lw=1.8, color='#555',
               label=f"TI 200 ps = {r['dG']:+.1f} kcal/mol")
ax_run.axhline(SHORT_TI[sol_main], ls=':', lw=1.8, color='#999',
               label=f"TI 40 ps  = {SHORT_TI[sol_main]:+.1f} kcal/mol")
ax_run.set_xlabel('λ', fontsize=13)
ax_run.set_ylabel('Cumulative ΔG (kcal/mol)', fontsize=13)
ax_run.set_title(f'(e) Running TI — {SOLVENT_LABELS[sol_main]}\nConvergence: 40 ps vs 200 ps',
                 fontsize=14, pad=8)
ax_run.legend(fontsize=11)
ax_run.grid(alpha=0.3)

# Row 1 right: OMP speedup ─────────────────────────────────────────────────
ax_omp = fig.add_subplot(gs[1, 2])
n_threads   = [1, 2, 4, 8]
sp_soft     = [1.0, 2.15, 4.04, 5.28]
sp_ideal    = [1.0, 2.0,  4.0,  8.0 ]
ax_omp.plot(n_threads, sp_ideal, 'k--', lw=1.8, label='Ideal linear', alpha=0.6)
ax_omp.plot(n_threads, sp_soft, 'o-', color='#E53935', lw=2.5, ms=9,
            label='lj/cut/soft/omp\n(32k atoms)')
ax_omp.fill_between(n_threads, sp_soft, 1, alpha=0.12, color='#E53935')

for nt, sp in zip(n_threads, sp_soft):
    ax_omp.annotate(f"{sp:.2f}×", (nt, sp), xytext=(nt + 0.15, sp - 0.45),
                    fontsize=12, color='#E53935', fontweight='bold')

ax_omp.set_xlabel('OMP threads', fontsize=13)
ax_omp.set_ylabel('Speedup vs. serial', fontsize=13)
ax_omp.set_title('(f) OpenMP FEP speedup\n(lj/cut/soft/omp, 32,000 atoms)',
                 fontsize=14, pad=8)
ax_omp.legend(fontsize=12, loc='upper left')
ax_omp.grid(alpha=0.3)
ax_omp.set_xticks(n_threads)
ax_omp.set_ylim(0.5, 9)

fig.suptitle(
    'Figure 6 · FEP Case Study: Li⁺ Desolvation in Battery Electrolyte Solvents\n'
    '(OPLS-AA force field; TI with 21 λ-windows; 40 ps and 200 ps/window)',
    fontsize=15, fontweight='bold', y=1.02)

out_pdf = str(OUTDIR / "figure6_fep_casestudy.pdf")
out_png = str(OUTDIR / "figure6_fep_casestudy.png")
plt.savefig(out_pdf, bbox_inches='tight', dpi=300)
plt.savefig(out_png, bbox_inches='tight', dpi=300)
print(f"Saved {out_pdf}")
print(f"       {out_png}")
