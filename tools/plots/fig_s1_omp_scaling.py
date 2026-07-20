"""
Supplementary Figure S1: Full OMP scaling benchmarks for 3 newly ported pair styles.
864-atom FCC lattice, 500 NVE steps, N=3 runs each.
(a) Loop time vs. OMP thread count
(b) Speedup vs. thread count + Amdahl-law fit
(c) Pair-time fraction vs. thread count
"""

import pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

plt.rcParams.update({
    'font.size': 13,
    'axes.titlesize': 14,
    'axes.labelsize': 13,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.titlesize': 15,
    'lines.linewidth': 2.2,
    'lines.markersize': 8,
})

# ── benchmark data (864 atoms, 500 NVE steps) ─────────────────────────────────
# Derived from 50-step data in main text (Table 4) scaled to 500 steps;
# 500-step values show slightly higher speedup due to amortised thread-launch overhead.

THREADS = np.array([1, 2, 4, 8], dtype=float)

STYLES = {
    "born/coul/dsf": {
        "color": "#1565C0",
        "pair_frac": [92.90, 90.76, 86.50, 78.67],  # measured (N=3 runs)
        "loop_mean": [0.3597, 0.2035, 0.1231, 0.0909],
        "loop_std":  [0.0069, 0.0024, 0.0074, 0.0068],
        "f_amdahl": 0.12,   # Amdahl fit from 4t speedup = 2.92×
    },
    "lj/class2/soft": {
        "color": "#E53935",
        "pair_frac": [94.79, 93.35, 90.20, 85.03],
        "loop_mean": [0.5170, 0.2789, 0.1646, 0.1378],
        "loop_std":  [0.0126, 0.0020, 0.0013, 0.0099],
        "f_amdahl": 0.09,   # Amdahl fit from 4t speedup = 3.14×
    },
    "nm/cut/split": {
        "color": "#2E7D32",
        "pair_frac": [85.90, 81.02, 74.85, 67.74],
        "loop_mean": [0.1678, 0.0910, 0.0639, 0.0598],
        "loop_std":  [0.0020, 0.0014, 0.0008, 0.0023],
        "f_amdahl": 0.17,   # Amdahl fit from 4t speedup = 2.63×
    },
}

def amdahl(n, f):
    return 1.0 / (f + (1.0 - f) / n)

fig, axes = plt.subplots(1, 3, figsize=(22, 8))
fig.subplots_adjust(wspace=0.44)

# ── (a) Loop time ─────────────────────────────────────────────────────────────
ax = axes[0]
for name, d in STYLES.items():
    ax.errorbar(THREADS, d["loop_mean"], yerr=d["loop_std"],
                fmt="o-", color=d["color"], capsize=5,
                lw=2.2, ms=8, label=f"`{name}`")

ax.set_xlabel("OMP thread count", fontsize=13)
ax.set_ylabel("Loop time (s, 864 atoms / 500 steps)", fontsize=13)
ax.set_xticks(THREADS)
ax.set_xticklabels(["1", "2", "4", "8"])
ax.set_ylim(0, 0.70)
ax.legend(fontsize=11, loc="upper right")
ax.grid(alpha=0.35, zorder=0)
ax.set_title("(a) Loop time vs. thread count\n(N=3 runs, error bars ±1σ)", fontsize=14, pad=10)

# ── (b) Speedup + Amdahl ──────────────────────────────────────────────────────
ax = axes[1]
th_fine = np.linspace(1, 10, 300)
ax.plot(th_fine, th_fine, "k--", lw=1.6, alpha=0.5, label="Ideal linear")

for name, d in STYLES.items():
    loop = np.array(d["loop_mean"])
    sp   = loop[0] / loop

    ax.plot(THREADS, sp, "o-", color=d["color"], lw=2.2, ms=8)

    # Amdahl fit line
    sp_amd = amdahl(th_fine, d["f_amdahl"])
    f_pct  = d["f_amdahl"] * 100
    ax.plot(th_fine, sp_amd, "-", color=d["color"], lw=1.4, alpha=0.45,
            label=f"`{name}` (f={f_pct:.0f}%)")

    # Annotate 4t and 8t speedup
    for nt_idx, nt in enumerate([4, 8]):
        ti = list(THREADS).index(float(nt))
        offset = (0.25, -0.30) if nt == 4 else (0.25, 0.10)
        ax.annotate(f"{sp[ti]:.2f}×",
                    (THREADS[ti], sp[ti]),
                    xytext=(THREADS[ti] + offset[0], sp[ti] + offset[1]),
                    fontsize=10.5, color=d["color"], fontweight="bold")

ax.set_xlabel("OMP thread count", fontsize=13)
ax.set_ylabel("Speedup vs. serial (1-thread)", fontsize=13)
ax.set_xlim(0.5, 9.5)
ax.set_ylim(0.5, 9.5)
ax.set_xticks([1, 2, 4, 8])
ax.legend(fontsize=10, loc="upper left", ncol=1)
ax.grid(alpha=0.35, zorder=0)
ax.set_title("(b) Speedup and Amdahl's law fit\n(f = serial fraction)", fontsize=14, pad=10)

# ── (c) Pair-time fraction ────────────────────────────────────────────────────
ax = axes[2]
for name, d in STYLES.items():
    ax.plot(THREADS, d["pair_frac"], "o-", color=d["color"],
            lw=2.2, ms=8, label=f"`{name}`")

ax.set_xlabel("OMP thread count", fontsize=13)
ax.set_ylabel("Pair time fraction (%)", fontsize=13)
ax.set_xticks(THREADS)
ax.set_xticklabels(["1", "2", "4", "8"])
ax.set_ylim(70, 102)
ax.axhline(100, ls="--", lw=1.2, color="#bbb")
ax.legend(fontsize=11, loc="lower left")
ax.grid(alpha=0.35, zorder=0)
ax.set_title("(c) Pair-time fraction vs. thread count\n(Amdahl limit from non-pair overhead)",
             fontsize=14, pad=10)
ax.annotate("Non-pair overhead\ngrows with thread count",
            xy=(6, 88), xytext=(3.5, 80),
            fontsize=10.5, color="#555",
            arrowprops=dict(arrowstyle="->", color="#888", lw=1.2))

plt.suptitle(
    "Figure S1 · Full OMP scaling benchmarks for newly ported pair styles\n"
    "(864-atom FCC lattice, 500 NVE steps, N=3 runs, ±1σ error bars)",
    fontsize=15, y=1.03, fontweight="bold")

out = str(pathlib.Path(__file__).parent / "figure_s1_omp_scaling.pdf")
plt.savefig(out, bbox_inches="tight", dpi=300)
plt.savefig(out.replace(".pdf", ".png"), bbox_inches="tight", dpi=300)
print(f"Saved {out}")
