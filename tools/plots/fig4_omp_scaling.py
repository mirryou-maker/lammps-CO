"""
Figure 4: OpenMP threading performance.
(a) Loop time vs thread count — 32k-atom lj/cut benchmark
(b) Speedup vs thread count with Amdahl's law fit
(c) Section fractions (Pair/Neigh/Comm/Other) vs thread count
(d) Small-system (864 atoms) speedup for 3 newly parallelized styles
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

# ── data ───────────────────────────────────────────────────────────────────
threads_large = np.array([1,     1,      4,     8    ])
labels_large  = ["Serial", "OMP 1t", "OMP 4t", "OMP 8t"]
loop_large    = np.array([4.333, 4.203, 1.004, 0.878])
std_large     = np.array([0.112, 0.102, 0.024, 0.117])
speedup_large = loop_large[0] / loop_large

pair_frac  = np.array([80.8, 80.5, 72.7, 69.0])
neigh_frac = np.array([15.7, 16.1, 16.8, 14.3])
comm_frac  = np.array([1.1,  1.1,  5.1,  7.0 ])
other_frac = 100 - pair_frac - neigh_frac - comm_frac

def amdahl_speedup(n, f):
    return 1.0 / (f + (1 - f) / n)

th_data = np.array([1, 4, 8], dtype=float)
sp_data = np.array([1.00, 4.31, 4.93])
popt, _ = curve_fit(amdahl_speedup, th_data, sp_data, p0=[0.18], bounds=(0, 1))
f_serial = popt[0]
th_fine  = np.linspace(1, 16, 200)
sp_amdahl = amdahl_speedup(th_fine, f_serial)
sp_ideal  = th_fine

new_styles  = ["nm/cut/split", "born/coul/dsf", "lj/class2/soft"]
serial_time = np.array([0.0332, 0.0463, 0.0934])
omp4_time   = np.array([0.0112, 0.0132, 0.0300])
speedup_new = serial_time / omp4_time

C_SERIAL = "#B0BEC5"; C_1T = "#90A4AE"; C_4T = "#1565C0"; C_8T = "#0D47A1"
bar_colors = [C_SERIAL, C_1T, C_4T, C_8T]

fig, axes = plt.subplots(2, 2, figsize=(20, 15))
fig.subplots_adjust(hspace=0.50, wspace=0.46)

# ── (a) loop time ──────────────────────────────────────────────────────────
ax = axes[0, 0]
x = np.arange(4)
bars = ax.bar(x, loop_large, yerr=std_large, capsize=6,
              color=bar_colors, edgecolor="white",
              error_kw={"elinewidth": 1.8, "ecolor": "#444"}, zorder=3)
for bar, val, sp in zip(bars, loop_large, speedup_large):
    idx = list(loop_large).index(val)
    ax.text(bar.get_x() + bar.get_width() / 2,
            val + std_large[idx] + 0.10,
            f"{val:.3f}s\n({sp:.2f}×)", ha="center", va="bottom",
            fontsize=11, color="#222", fontweight="bold")
ax.set_xticks(x)
ax.set_xticklabels(labels_large, fontsize=12)
ax.set_ylabel("Loop time (s)", fontsize=13)
ax.set_ylim(0, 5.8)
ax.set_title("(a) Loop time vs. thread count\n(32,000 atoms, 500 steps, lj/cut)", fontsize=14, pad=10)
ax.grid(axis="y", alpha=0.35, zorder=0)

# ── (b) speedup with Amdahl ────────────────────────────────────────────────
ax = axes[0, 1]
ax.plot(th_fine, sp_ideal,   "--", color="#bbb",    lw=1.8, label="Ideal (linear)")
ax.plot(th_fine, sp_amdahl,  "-",  color="#E53935", lw=2.2,
        label=f"Amdahl fit (serial frac = {f_serial:.2f})")
ax.scatter([1, 4, 8], [1.00, 4.31, 4.93], s=100, zorder=5, color=C_4T,
           label="Measured", edgecolors="white", linewidths=1.5)

offsets = [(0.5, -0.35), (0.5, -0.35), (0.5, -0.35)]
for th, sp, lab, (dx, dy) in zip([1, 4, 8], [1.00, 4.31, 4.93],
                                  ["1.00×", "4.31×", "4.93×"], offsets):
    ax.annotate(lab, (th, sp), xytext=(th + dx, sp + dy),
                fontsize=12, color=C_4T, fontweight="bold")

ax.set_xlabel("OMP thread count", fontsize=13)
ax.set_ylabel("Speedup vs. serial baseline", fontsize=13)
ax.set_xlim(0.5, 10.5)
ax.set_ylim(0.5, 10.5)
ax.set_xticks([1, 2, 4, 8, 10])
ax.legend(fontsize=11, loc="upper left")
ax.set_title("(b) Speedup scaling and Amdahl's law fit", fontsize=14, pad=10)
ax.grid(alpha=0.35, zorder=0)

# ── (c) section fractions ──────────────────────────────────────────────────
ax = axes[1, 0]
x = np.arange(4)
w = 0.6
ax.bar(x, pair_frac,  w, label="Pair",  color="#1565C0", zorder=3)
ax.bar(x, neigh_frac, w, bottom=pair_frac, label="Neigh", color="#43A047", zorder=3)
ax.bar(x, comm_frac,  w, bottom=pair_frac + neigh_frac, label="Comm", color="#FB8C00", zorder=3)
ax.bar(x, other_frac, w, bottom=pair_frac + neigh_frac + comm_frac,
       label="Other", color="#90A4AE", zorder=3)
ax.set_xticks(x)
ax.set_xticklabels(labels_large, fontsize=12)
ax.set_ylabel("Fraction of loop time (%)", fontsize=13)
ax.set_ylim(0, 106)
ax.legend(fontsize=12, loc="lower right")
ax.set_title("(c) Section fractions vs. thread count\n(Amdahl limit from non-pair sections)",
             fontsize=14, pad=10)
ax.grid(axis="y", alpha=0.35, zorder=0)

# ── (d) new-style 4-thread speedup ────────────────────────────────────────
ax = axes[1, 1]
x2 = np.arange(3)
w2 = 0.35
bser = ax.bar(x2 - w2 / 2, serial_time * 1000, w2, color="#B0BEC5", label="Serial (1t)", zorder=3)
bomp = ax.bar(x2 + w2 / 2, omp4_time   * 1000, w2, color="#1565C0", label="OMP 4t",     zorder=3)

for bar, sp in zip(bomp, speedup_new):
    ax.text(bar.get_x() + bar.get_width() / 2 + 0.02,
            bar.get_height() + 1.0,
            f"{sp:.2f}×", ha="center", va="bottom", fontsize=13,
            color="#1565C0", fontweight="bold")

ax.set_xticks(x2)
ax.set_xticklabels(new_styles, fontsize=12)
ax.set_ylabel("Loop time (ms, 864 atoms / 50 steps)", fontsize=12)
ax.legend(fontsize=12)
ax.set_title("(d) Speedup for 3 newly parallelized styles\n(4 OMP threads vs. serial)",
             fontsize=14, pad=10)
ax.grid(axis="y", alpha=0.35, zorder=0)

plt.suptitle("Figure 4 · OpenMP performance benchmarks", fontsize=15, y=1.02, fontweight="bold")

out = str(pathlib.Path(__file__).parent / "figure4_omp_scaling.pdf")
plt.savefig(out, bbox_inches="tight", dpi=300)
plt.savefig(out.replace(".pdf", ".png"), bbox_inches="tight", dpi=300)
print(f"Saved {out}")
