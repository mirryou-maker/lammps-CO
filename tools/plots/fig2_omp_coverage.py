"""
Figure 2: OpenMP acceleration coverage in LAMMPS before and after optimization.
(a) Bar chart per interaction category
(b) Pie chart of 35 new variants by source package
(c) Radar chart comparing OMP / INTEL / OPT / KOKKOS coverage
"""

import pathlib
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
    'lines.linewidth': 2.2,
})

C_BEFORE = "#B0BEC5"
C_AFTER  = "#1565C0"
C_EXTRA  = "#1E88E5"
C_DIEL   = "#43A047"
C_KSPACE = "#FB8C00"

categories  = ["pair\n(307)", "bond\n(28)", "angle\n(30)", "dihedral\n(21)", "improper\n(15)"]
total       = [307, 28, 30, 21, 15]
before      = [157, 16, 23, 17, 11]
after       = [192, 16, 23, 17, 11]
before_pct  = [100 * b / t for b, t in zip(before, total)]
after_pct   = [100 * a / t for a, t in zip(after,  total)]

x = np.arange(len(categories))
w = 0.35

fig, axes = plt.subplots(1, 3, figsize=(22, 8))
fig.subplots_adjust(wspace=0.52)

# ── (a) coverage per category ──────────────────────────────────────────────
ax = axes[0]
bars_b = ax.bar(x - w / 2, before_pct, w, color=C_BEFORE, label="Before", zorder=3)
bars_a = ax.bar(x + w / 2, after_pct,  w, color=C_AFTER,  label="After",  zorder=3)

for bar, pct in zip(bars_b, before_pct):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1.5,
            f"{pct:.0f}%", ha="center", va="bottom", fontsize=11, color="#555")
for bar, pct in zip(bars_a, after_pct):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 5.5,
            f"{pct:.0f}%", ha="center", va="bottom", fontsize=11,
            color=C_AFTER, fontweight="bold")

ax.set_xticks(x)
ax.set_xticklabels(categories, fontsize=12)
ax.set_ylabel("OpenMP coverage (%)", fontsize=13)
ax.set_ylim(0, 110)
ax.axhline(100, ls="--", lw=1.2, color="#aaa")
ax.legend(fontsize=12, loc="lower right")
ax.grid(axis="y", alpha=0.35, zorder=0)
ax.set_title("(a) Coverage by interaction category", fontsize=14, pad=10)

# ── (b) pie chart of 35 new styles by package ─────────────────────────────
ax = axes[1]
labels  = ["EXTRA-PAIR\n(18 styles)", "DIELECTRIC\n(7 styles)",
           "KSPACE\n(3 styles)",      "Other EXTRA\n(7 styles)"]
sizes   = [18, 7, 3, 7]
colors  = [C_EXTRA, C_DIEL, C_KSPACE, "#AB47BC"]
explode = [0.04] * 4

wedges, texts, autotexts = ax.pie(
    sizes, explode=explode, labels=labels, colors=colors,
    autopct="%1.0f%%", startangle=105,
    textprops={"fontsize": 11},
    wedgeprops={"linewidth": 1.0, "edgecolor": "white"},
    pctdistance=0.72, labeldistance=1.18)
for at in autotexts:
    at.set_fontsize(12)
    at.set_fontweight("bold")
ax.set_title("(b) 35 new OMP variants by package", fontsize=14, pad=10)

# ── (c) radar chart (OMP / INTEL / OPT / KOKKOS) ─────────────────────────
ax = axes[2]
ax.remove()
ax = fig.add_subplot(1, 3, 3, polar=True)

radar_cats = ["pair", "bond", "angle", "dihedral", "improper"]
n_cats     = len(radar_cats)
angles     = np.linspace(0, 2 * np.pi, n_cats, endpoint=False).tolist()
angles    += angles[:1]

packages = {
    "OMP (after)"  : ([63, 57, 77, 81, 73], C_AFTER,   "-",  2.8),
    "OMP (before)" : ([51, 57, 77, 81, 73], C_BEFORE,  "--", 1.8),
    "INTEL"        : ([6,  0,  0,  0,  0 ], "#FF7043",  "-",  1.8),
    "OPT"          : ([5,  0,  0,  0,  0 ], "#66BB6A",  "-",  1.8),
    "KOKKOS"       : ([33, 50, 60, 76, 67], "#FFA726",  "-",  1.8),
}

for name, (vals, col, ls, lw) in packages.items():
    v = vals + vals[:1]
    ax.plot(angles, v, color=col, linewidth=lw, linestyle=ls, label=name)
    ax.fill(angles, v, alpha=0.07, color=col)

ax.set_xticks(angles[:-1])
ax.set_xticklabels(radar_cats, fontsize=12)
ax.set_yticks([25, 50, 75, 100])
ax.set_yticklabels(["25%", "50%", "75%", "100%"], fontsize=10, color="#777")
ax.set_ylim(0, 100)
ax.legend(loc="upper right", bbox_to_anchor=(1.55, 1.18), fontsize=11)
ax.set_title("(c) Coverage across\nacceleration packages", fontsize=14, pad=18)

plt.suptitle(
    "Figure 2 · OpenMP pair-style coverage before and after AI-assisted optimization",
    fontsize=15, y=1.03, fontweight="bold")

out = str(pathlib.Path(__file__).parent / "figure2_omp_coverage.pdf")
plt.savefig(out, bbox_inches="tight", dpi=300)
plt.savefig(out.replace(".pdf", ".png"), bbox_inches="tight", dpi=300)
print(f"Saved {out}")
