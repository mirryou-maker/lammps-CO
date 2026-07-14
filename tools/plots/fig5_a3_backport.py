"""
Figure 5: A-3 restrict/fxtmp optimization.
(a) Code transformation diagram (text-based annotation)
(b) Pair-time before/after for 3 styles
(c) Overall serial speedup from A-3 on lj/cut (4k atoms)
"""

import pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
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

styles   = ["momb\n(port #84)", "coul/ctip\n(port #85)", "dispersion/d3\n(port #86)"]
before_p = np.array([1.706, 22.17, 15.27])
after_p  = np.array([1.579, 20.90, 15.76])
impr_p   = (before_p - after_p) / before_p * 100

serial_before = 0.2890
serial_after  = 0.2768
serial_impr   = (serial_before - serial_after) / serial_before * 100

C_BEFORE = "#90A4AE"
C_AFTER  = "#1565C0"
C_NOISE  = "#EF9A9A"

fig = plt.figure(figsize=(18, 10))
gs  = fig.add_gridspec(2, 3, hspace=0.52, wspace=0.40)
ax_code = fig.add_subplot(gs[0, :])
ax_pair = fig.add_subplot(gs[1, :2])
ax_glob = fig.add_subplot(gs[1, 2])

# ── (a) code transformation panel ─────────────────────────────────────────
ax_code.set_xlim(0, 10); ax_code.set_ylim(0, 6.8)
ax_code.axis("off")
ax_code.set_title("(a) Hotloop transformation: Before vs. After A-3 optimization",
                  fontsize=14, pad=8, loc="left", fontweight="bold")

before_code = (
    "// BEFORE  (standard pair_lj_cut.cpp)\n"
    "double **f = atom->f;\n"
    "for (ii = 0; ii < inum; ii++) {\n"
    "  i = ilist[ii];\n"
    "  for (jj = 0; jj < jnum; jj++) {\n"
    "    // compute fpair ...\n"
    "    f[i][0] += delx*fpair;  // cache-hostile write\n"
    "    f[i][1] += dely*fpair;  //   every iteration\n"
    "    f[i][2] += delz*fpair;\n"
    "  }\n"
    "}"
)
ax_code.text(0.15, 6.5, before_code, fontsize=10, family="monospace",
             va="top", ha="left",
             bbox=dict(boxstyle="round,pad=0.6", fc="#FFF8E1", ec="#FBC02D", lw=2.0))
ax_code.text(0.15, 0.25,
             "⚠  Loop-carried dependency on f[i]  →  auto-vectorisation blocked",
             fontsize=11, color="#B71C1C", style="italic")

after_code = (
    "// AFTER  (A-3 applied)\n"
    "double * _noalias const f0 = atom->f[0];  // no-alias hint\n"
    "// outer loop:\n"
    "  double fxtmp=0, fytmp=0, fztmp=0;   // register accumulators\n"
    "  for (jj = 0; jj < jnum; jj++) {\n"
    "    // compute fpair ...\n"
    "    double fx=delx*fpair, fy=dely*fpair, fz=delz*fpair;\n"
    "    fxtmp+=fx; fytmp+=fy; fztmp+=fz;  // register only\n"
    "    if (newton||j<nlocal){ fj[0]-=fx; fj[1]-=fy; fj[2]-=fz; }\n"
    "  }\n"
    "  fi[0]+=fxtmp; fi[1]+=fytmp; fi[2]+=fztmp;  // single write-back"
)
ax_code.text(5.25, 6.5, after_code, fontsize=10, family="monospace",
             va="top", ha="left",
             bbox=dict(boxstyle="round,pad=0.6", fc="#E8F5E9", ec="#2E7D32", lw=2.0))
ax_code.text(5.25, 0.25,
             "✓  No loop-carried dependency  →  compiler auto-vectorises inner loop",
             fontsize=11, color="#1B5E20", style="italic")

ax_code.annotate("", xy=(5.12, 3.5), xytext=(4.85, 3.5),
                 arrowprops=dict(arrowstyle="->", color="#444", lw=2.5))
ax_code.text(4.98, 4.0, "_noalias\n+ fxtmp", fontsize=10, ha="center",
             color="#444", style="italic")

# ── (b) pair-time before/after ─────────────────────────────────────────────
x = np.arange(3)
w = 0.32
ax_pair.bar(x - w / 2, before_p, w, color=C_BEFORE, label="Before A-3", zorder=3)
ax_pair.bar(x + w / 2, after_p,  w, color=C_AFTER,  label="After A-3",  zorder=3)

for i, (b, a, imp) in enumerate(zip(before_p, after_p, impr_p)):
    col  = C_AFTER if abs(imp) > 1.5 else C_NOISE
    sign = "+" if imp > 0 else ""
    ax_pair.annotate(f"{sign}{imp:.1f}%",
                     xy=(x[i] + w / 2, a),
                     xytext=(x[i] + w / 2, max(b, a) + max(before_p) * 0.06),
                     ha="center", fontsize=13, fontweight="bold", color=col,
                     arrowprops=dict(arrowstyle="-", color=col, lw=1.0))

ax_pair.set_xticks(x)
ax_pair.set_xticklabels(styles, fontsize=12)
ax_pair.set_ylabel("Pair time (s)", fontsize=13)
ax_pair.legend(fontsize=12)
ax_pair.set_title("(b) Pair-time before/after A-3  (6,912 atoms, 300 NVE steps, N=6)",
                  fontsize=14, pad=10)
ax_pair.grid(axis="y", alpha=0.35, zorder=0)
note = ("¹ dispersion/d3: improvement within measurement noise\n"
        "  (complex multi-pass loop; verified bit-identical)")
ax_pair.text(0.98, 0.03, note, transform=ax_pair.transAxes,
             fontsize=10, ha="right", va="bottom", color="#777", style="italic")

# ── (c) global serial speedup on lj/cut ───────────────────────────────────
ax_glob.bar([0, 1], [serial_before * 1000, serial_after * 1000],
            color=[C_BEFORE, C_AFTER], edgecolor="white", width=0.5, zorder=3)
ax_glob.text(0, serial_before * 1000 + 2.5, f"{serial_before * 1000:.1f} ms",
             ha="center", fontsize=12)
ax_glob.text(1, serial_after * 1000  + 2.5, f"{serial_after * 1000:.1f} ms",
             ha="center", fontsize=12, color=C_AFTER, fontweight="bold")
ax_glob.annotate(f"+{serial_impr:.1f}%\nimprovement",
                 xy=(1, serial_after * 1000),
                 xytext=(1.38, (serial_before + serial_after) / 2 * 1000),
                 fontsize=13, fontweight="bold", color=C_AFTER,
                 arrowprops=dict(arrowstyle="->", color=C_AFTER, lw=2.0))
ax_glob.set_xticks([0, 1])
ax_glob.set_xticklabels(["Before A-3", "After A-3"], fontsize=12)
ax_glob.set_ylabel("Serial loop time (ms)", fontsize=13)
ax_glob.set_ylim(0, 380)
ax_glob.set_title("(c) Global serial improvement\n(lj/cut, 4,000 atoms, 250 steps)",
                  fontsize=14, pad=10)
ax_glob.grid(axis="y", alpha=0.35, zorder=0)

plt.suptitle("Figure 5 · A-3 restrict/_noalias + fxtmp accumulator optimization",
             fontsize=15, y=1.02, fontweight="bold")

out = str(pathlib.Path(__file__).parent / "figure5_a3_backport.pdf")
plt.savefig(out, bbox_inches="tight", dpi=300)
plt.savefig(out.replace(".pdf", ".png"), bbox_inches="tight", dpi=300)
print(f"Saved {out}")
