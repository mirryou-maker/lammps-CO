"""
Figure 1: AI-assisted LAMMPS optimization workflow diagram.
(a) Five-step pipeline with risk/LLM-suitability axes
(b) Claude Code interaction pattern
(c) Three-level verification pyramid
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
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

RISK_COLS = ["#43A047", "#66BB6A", "#1E88E5", "#1565C0", "#0D47A1"]

fig, axes = plt.subplots(1, 3, figsize=(20, 7))
fig.subplots_adjust(wspace=0.18)

# ── (a) five-step pipeline ─────────────────────────────────────────────────
ax = axes[0]
ax.set_xlim(0, 10); ax.set_ylim(-0.5, 12); ax.axis("off")
ax.set_title("(a) Five-step optimization pipeline", fontsize=14, pad=10, fontweight="bold")

steps = [
    ("B-1", "Benchmarking\nInfrastructure",  "No code changes\nAutomates measurement",   10.0, "Very Low"),
    ("B-2", "Coverage\nDiagnostics",          "Script analysis\nno source changes",          8.1, "Very Low"),
    ("A-2", "Compiler Flag\nTuning",           "CMake preset\nno source changes",             6.2, "Low"),
    ("A-1", "OMP Variant\nGeneration",         "35 new files\nbit-identical verified",        4.3, "Low"),
    ("A-3", "restrict/fxtmp\nBackport",        "~102 files modified\nbit-identical OK",       2.4, "Low–Med"),
]

for i, (tag, title, desc, y, risk) in enumerate(steps):
    col = RISK_COLS[i]
    rect = FancyBboxPatch((0.4, y - 0.90), 9.0, 1.65,
                          boxstyle="round,pad=0.18", linewidth=2.0,
                          edgecolor=col, facecolor=col + "22")
    ax.add_patch(rect)
    ax.text(1.1, y - 0.05, tag, fontsize=13, fontweight="bold", color=col, va="center")
    ax.text(3.1, y + 0.30, title, fontsize=12, fontweight="bold", color="#222", va="center")
    ax.text(3.1, y - 0.45, desc, fontsize=10.5, color="#444", va="center")
    ax.text(9.3, y - 0.05, f"Risk: {risk}", fontsize=10, color=col, va="center", ha="right")
    if i < len(steps) - 1:
        ax.annotate("", xy=(5, steps[i + 1][3] + 0.75), xytext=(5, y - 0.90),
                    arrowprops=dict(arrowstyle="-|>", color="#888", lw=1.5))

ax.text(5, 11.6, "Increasing risk & code invasiveness →",
        ha="center", fontsize=10.5, color="#777", style="italic")

# ── (b) Claude Code interaction pattern ───────────────────────────────────
ax = axes[1]
ax.set_xlim(0, 10); ax.set_ylim(0, 10.5); ax.axis("off")
ax.set_title("(b) Claude Code interaction pattern", fontsize=14, pad=10, fontweight="bold")

boxes = [
    (5.0, 9.3, "User specifies\noptimization target",   "#E3F2FD", "#1565C0"),
    (5.0, 7.4, "Claude Code reads\nsource + reference", "#F3E5F5", "#7B1FA2"),
    (5.0, 5.5, "Pattern extraction\n& code generation",  "#E8F5E9", "#2E7D32"),
    (5.0, 3.6, "Automated\nnumerical verification",      "#FFF3E0", "#E65100"),
    (5.0, 1.7, "Pass → commit\nFail → revise prompt",   "#E8EAF6", "#283593"),
]

for (x, y, txt, fc, ec) in boxes:
    rect = FancyBboxPatch((x - 3.9, y - 0.75), 7.8, 1.38,
                          boxstyle="round,pad=0.14", linewidth=2.0,
                          edgecolor=ec, facecolor=fc)
    ax.add_patch(rect)
    ax.text(x, y - 0.05, txt, ha="center", va="center", fontsize=11.5, color="#222")

for i in range(len(boxes) - 1):
    ax.annotate("", xy=(5, boxes[i + 1][1] + 0.63), xytext=(5, boxes[i][1] - 0.75),
                arrowprops=dict(arrowstyle="-|>", color="#666", lw=1.6))

ax.annotate("", xy=(1.0, boxes[2][1] + 0.20), xytext=(1.0, boxes[4][1]),
            arrowprops=dict(arrowstyle="-|>", color="#E65100", lw=1.8,
                            connectionstyle="arc3,rad=-0.35"))
ax.text(0.0, 3.2, "Revise\nprompt", fontsize=10, color="#E65100",
        ha="center", va="center", rotation=90)

# ── (c) verification pyramid ──────────────────────────────────────────────
ax = axes[2]
ax.set_xlim(0, 10); ax.set_ylim(-0.5, 10.5); ax.axis("off")
ax.set_title("(c) Three-level verification pyramid", fontsize=14, pad=10, fontweight="bold")

levels = [
    (0.4,  2.4, 4.8, "Level 1: Build-time",
     "Compilation succeeds\nno new MSVC warnings", "#66BB6A"),
    (3.0,  2.4, 3.4, "Level 2: Bit-identical",
     "50 NVE steps: temp, epair,\netotal, press match to 15 sf", "#1E88E5"),
    (5.6,  2.4, 2.1, "Level 3: Conservation",
     "250 NVE steps\nenergy drift < 0.01%", "#1565C0"),
]

for (yb, h, wh, lab, sub, col) in levels:
    tx = [5 - wh, 5 + wh, 5, 5 - wh]
    ty = [yb, yb, yb + h, yb]
    ax.fill(tx, ty, color=col + "33", zorder=2)
    ax.plot(tx, ty, color=col, lw=2.2, zorder=3)
    ax.text(5, yb + h * 0.40, lab, ha="center", va="center",
            fontsize=11, fontweight="bold", color=col)
    ax.text(5, yb + h * 0.10, sub, ha="center", va="center",
            fontsize=9.5, color="#333")

scope_labels = ["All 137 files", "35 OMP variants", "8 styles"]
for (yb, h, _, _, _, _), n in zip(levels, scope_labels):
    ax.text(9.6, yb + h / 2, n, ha="right", va="center",
            fontsize=10, color="#555", style="italic")
ax.text(9.7, 8.0, "Scope", ha="right", fontsize=11, color="#777", fontweight="bold")

ax.text(5, 8.2, "✓ All 35 new OMP variants\npass all 3 levels",
        ha="center", va="center", fontsize=11, color="#0D47A1",
        bbox=dict(fc="#E8EAF6", ec="#3949AB", boxstyle="round,pad=0.5", lw=1.8))

plt.suptitle(
    "Figure 1 · AI-assisted LAMMPS optimization: workflow, interaction pattern, and verification",
    fontsize=15, y=1.02, fontweight="bold")

out = "D:/Claude-Code-R/LAMMPS-CO/tools/plots/figure1_workflow.pdf"
plt.savefig(out, bbox_inches="tight", dpi=300)
plt.savefig(out.replace(".pdf", ".png"), bbox_inches="tight", dpi=300)
print(f"Saved {out}")
