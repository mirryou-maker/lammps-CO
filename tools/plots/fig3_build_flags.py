"""
Figure 3: Compiler flag optimization results.
(a) Loop time (mean ± 1σ) for 5 build configs — large system (32k atoms)
(b) % improvement vs baseline
(c) Pair time fraction stays constant → flags don't hurt non-pair sections
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
    'lines.linewidth': 2.0,
})

configs   = ["Baseline\n(no flags)", "/arch:AVX2", "/fp:fast", "LTO\n(/GL+LTCG)", "AVX2 +\n/fp:fast"]
loop_mean = np.array([4.333, 4.209, 4.214, 4.245, 4.058])
loop_std  = np.array([0.112, 0.169, 0.043, 0.199, 0.124])
pair_pct  = np.array([80.8, 80.5, 80.3, 80.5, 80.3])
pct_impr  = (loop_mean[0] - loop_mean) / loop_mean[0] * 100

colors = ["#90A4AE", "#1565C0", "#1E88E5", "#E53935", "#0D47A1"]

fig, axes = plt.subplots(1, 3, figsize=(20, 8))
fig.subplots_adjust(wspace=0.48)

x = np.arange(len(configs))

# ── (a) loop time ──────────────────────────────────────────────────────────
ax = axes[0]
bars = ax.bar(x, loop_mean, yerr=loop_std, capsize=6,
              color=colors, edgecolor="white", linewidth=0.9,
              error_kw={"elinewidth": 1.8, "ecolor": "#444"}, zorder=3)
ax.axhline(loop_mean[0], ls="--", lw=1.4, color="#aaa", label="Baseline")
for bar, val in zip(bars, loop_mean):
    ax.text(bar.get_x() + bar.get_width() / 2, val + 0.15,
            f"{val:.3f}s", ha="center", va="bottom", fontsize=11)
ax.set_xticks(x)
ax.set_xticklabels(configs, fontsize=11.5)
ax.set_ylabel("Loop time (s)", fontsize=13)
ax.set_ylim(0, 5.5)
ax.set_title("(a) Loop time — 32,000-atom benchmark", fontsize=14, pad=10)
ax.grid(axis="y", alpha=0.35, zorder=0)

# ── (b) % improvement ─────────────────────────────────────────────────────
ax = axes[1]
bar_colors = []
for p in pct_impr:
    if p < -0.5:
        bar_colors.append("#E53935")
    elif p < 0.5:
        bar_colors.append("#90A4AE")
    else:
        bar_colors.append("#1565C0")

bars2 = ax.bar(x, pct_impr, color=bar_colors, edgecolor="white", linewidth=0.9, zorder=3)
ax.axhline(0, lw=1.2, color="#444")
for bar, val in zip(bars2, pct_impr):
    offset = 0.25 if val >= 0 else -0.60
    ax.text(bar.get_x() + bar.get_width() / 2, val + offset,
            f"{val:+.1f}%", ha="center", va="bottom", fontsize=11.5, fontweight="bold")
ax.set_xticks(x)
ax.set_xticklabels(configs, fontsize=11.5)
ax.set_ylabel("Loop-time improvement vs. baseline (%)", fontsize=12.5)
ax.set_ylim(-3.5, 11)
ax.set_title("(b) Percentage improvement", fontsize=14, pad=10)
ax.grid(axis="y", alpha=0.35, zorder=0)

pos_patch = mpatches.Patch(color="#1565C0", label="Speed-up")
neg_patch = mpatches.Patch(color="#E53935", label="Slow-down")
ax.legend(handles=[pos_patch, neg_patch], fontsize=12)

# ── (c) pair fraction ─────────────────────────────────────────────────────
ax = axes[2]
bars3 = ax.bar(x, pair_pct, color=colors, edgecolor="white", linewidth=0.9, zorder=3)
ax.axhline(pair_pct[0], ls="--", lw=1.4, color="#aaa")
for bar, val in zip(bars3, pair_pct):
    ax.text(bar.get_x() + bar.get_width() / 2, val + 0.10,
            f"{val:.1f}%", ha="center", va="bottom", fontsize=11)
ax.set_xticks(x)
ax.set_xticklabels(configs, fontsize=11.5)
ax.set_ylabel("Pair time fraction (%)", fontsize=13)
ax.set_ylim(78.5, 83.0)
ax.set_title("(c) Pair fraction — unaffected by flags", fontsize=14, pad=10)
ax.grid(axis="y", alpha=0.35, zorder=0)

plt.suptitle(
    "Figure 3 · Compiler flag optimization on 32,000-atom FCC Lennard-Jones benchmark (N=5 runs)",
    fontsize=15, y=1.02, fontweight="bold")

out = str(pathlib.Path(__file__).parent / "figure3_build_flags.pdf")
plt.savefig(out, bbox_inches="tight", dpi=300)
plt.savefig(out.replace(".pdf", ".png"), bbox_inches="tight", dpi=300)
print(f"Saved {out}")
