"""
Figure 6: Combined optimization summary.
(a) Waterfall chart: serial → AVX2 → A-3 → OMP 4t → OMP 8t
(b) Tornado chart: individual % contributions
(c) 2D heatmap: flags × OMP threads → speedup
(d) OMP coverage donut: before/after with 35-new highlight
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

plt.rcParams.update({
    'font.size': 13,
    'axes.titlesize': 14,
    'axes.labelsize': 13,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.titlesize': 15,
    'lines.linewidth': 2.0,
})

BASELINE = 4.333

steps = [
    ("Serial\nbaseline",  4.333, "#90A4AE"),
    ("+ AVX2\nflag",      4.209, "#1976D2"),
    ("+ A-3\npatch",      4.058, "#0288D1"),
    ("OMP 4t\n(-sf omp)", 1.004, "#1565C0"),
    ("OMP 8t\n(-sf omp)", 0.878, "#0D47A1"),
]

flags_labels   = ["No flags", "AVX2", "AVX2+fp:fast"]
threads_labels = ["Serial", "1t", "4t", "8t"]
heat_data = np.array([
    [4.333, 4.203, 1.004, 0.878],
    [4.209, 4.000, 1.012, 0.950],
    [4.058, 3.850, 1.015, 0.947],
])
heat_speedup = BASELINE / heat_data

fig, axes = plt.subplots(2, 2, figsize=(16, 13))
fig.subplots_adjust(hspace=0.50, wspace=0.46)

# ── (a) waterfall ──────────────────────────────────────────────────────────
ax = axes[0, 0]
vals   = [s[1] for s in steps]
colors = [s[2] for s in steps]
labels = [s[0] for s in steps]
speedups = [BASELINE / v for v in vals]

x = np.arange(len(steps))
bars = ax.bar(x, vals, color=colors, edgecolor="white", linewidth=1.0, width=0.6, zorder=3)

for i, (bar, val, sp) in enumerate(zip(bars, vals, speedups)):
    ax.text(bar.get_x() + bar.get_width() / 2, val + 0.10,
            f"{val:.3f}s\n({sp:.2f}×)", ha="center", va="bottom",
            fontsize=11, fontweight="bold", color=colors[i])
    if i > 0:
        delta = vals[i - 1] - val
        pct   = delta / vals[i - 1] * 100
        ax.annotate(f"−{pct:.1f}%",
                    xy=(x[i], val + (vals[i - 1] - val) / 2 + 0.05),
                    fontsize=10, ha="center", color="#555", style="italic")

ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=11.5)
ax.set_ylabel("Loop time (s)", fontsize=13)
ax.set_ylim(0, 5.8)
ax.set_title("(a) Cumulative optimization waterfall\n(32,000 atoms, 500 steps)",
             fontsize=14, pad=10)
ax.grid(axis="y", alpha=0.30, zorder=0)

# ── (b) tornado chart ──────────────────────────────────────────────────────
ax = axes[0, 1]
tornado_items = [
    ("A-3 serial patch\n(lj/cut, 4k atoms)",   4.2,  "#26C6DA"),
    ("/arch:AVX2\n(32k atoms)",                  2.9,  "#1E88E5"),
    ("AVX2 + /fp:fast\n(32k atoms)",             6.3,  "#0D47A1"),
    ("OMP 4 threads\n(-sf omp, 32k atoms)", (BASELINE - 1.004) / BASELINE * 100, "#43A047"),
    ("OMP 8 threads\n(-sf omp, 32k atoms)", (BASELINE - 0.878) / BASELINE * 100, "#1B5E20"),
]
tlabels = [t[0] for t in tornado_items]
tvals   = [t[1] for t in tornado_items]
tcolors = [t[2] for t in tornado_items]

y = np.arange(len(tornado_items))
hbars = ax.barh(y, tvals, color=tcolors, edgecolor="white", height=0.6, zorder=3)
for bar, val in zip(hbars, tvals):
    ax.text(val + 0.8, bar.get_y() + bar.get_height() / 2,
            f"+{val:.1f}%", va="center", fontsize=12, fontweight="bold")
ax.set_yticks(y)
ax.set_yticklabels(tlabels, fontsize=11.5)
ax.set_xlabel("Loop-time reduction vs. serial baseline (%)", fontsize=12.5)
ax.set_xlim(0, 90)
ax.axvline(0, lw=1.0, color="#444")
ax.set_title("(b) Individual optimization contributions\n(% loop-time reduction)",
             fontsize=14, pad=10)
ax.grid(axis="x", alpha=0.35, zorder=0)

# ── (c) heatmap ────────────────────────────────────────────────────────────
ax = axes[1, 0]
im = ax.imshow(heat_speedup, cmap="YlOrRd", vmin=1.0, vmax=5.5, aspect="auto")
cbar = plt.colorbar(im, ax=ax)
cbar.set_label("Speedup vs. serial no-flags baseline", fontsize=12)
cbar.ax.tick_params(labelsize=11)

for i in range(heat_speedup.shape[0]):
    for j in range(heat_speedup.shape[1]):
        val = heat_speedup[i, j]
        tc  = "white" if val > 3.5 else "black"
        ax.text(j, i, f"{val:.2f}×", ha="center", va="center",
                fontsize=13, color=tc, fontweight="bold")

ax.set_xticks(range(4))
ax.set_xticklabels(threads_labels, fontsize=12)
ax.set_yticks(range(3))
ax.set_yticklabels(flags_labels, fontsize=12)
ax.set_xlabel("OMP thread count", fontsize=13)
ax.set_ylabel("Build flags", fontsize=13)

best_per_col = np.argmax(heat_speedup, axis=0)
for j, i in enumerate(best_per_col):
    ax.add_patch(plt.Rectangle((j - 0.48, i - 0.48), 0.96, 0.96,
                                fill=False, edgecolor="#222", linewidth=2.8))

ax.set_title("(c) Speedup heatmap — flags × threads\n(■ = Pareto-optimal config)",
             fontsize=14, pad=10)

# ── (d) coverage donut ─────────────────────────────────────────────────────
ax = axes[1, 1]
before_n = 157
new_n    = 35
skip_n   = 307 - 157 - 35

sizes   = [before_n, new_n, skip_n]
colors2 = ["#90A4AE", "#E53935", "#CFD8DC"]
labels2 = [f"Pre-existing OMP\n({before_n} styles)",
           f"Newly added\n({new_n} styles, this work)",
           f"Skip/other\n({skip_n} styles)"]
explode = [0, 0.08, 0]

wedges, texts, autotexts = ax.pie(
    sizes, explode=explode, labels=None, colors=colors2,
    autopct="%1.0f%%", startangle=90,
    wedgeprops={"linewidth": 1.2, "edgecolor": "white"},
    pctdistance=0.72)
for at in autotexts:
    at.set_fontsize(13)
    at.set_fontweight("bold")

ax.text(0, 0, "307\npair\nstyles", ha="center", va="center",
        fontsize=12, fontweight="bold", color="#333")

patches = [mpatches.Patch(color=c, label=l) for c, l in zip(colors2, labels2)]
ax.legend(handles=patches, fontsize=11, loc="lower center",
          bbox_to_anchor=(0.5, -0.26), ncol=1)
ax.set_title("(d) pair_style OMP coverage distribution\n(51.1% → 62.5%, +11.4 pp)",
             fontsize=14, pad=10)

plt.suptitle("Figure 6 · Combined optimization summary", fontsize=15, y=1.02, fontweight="bold")

out = "D:/Claude-Code-R/LAMMPS-CO/tools/plots/figure6_summary.pdf"
plt.savefig(out, bbox_inches="tight", dpi=300)
plt.savefig(out.replace(".pdf", ".png"), bbox_inches="tight", dpi=300)
print(f"Saved {out}")
