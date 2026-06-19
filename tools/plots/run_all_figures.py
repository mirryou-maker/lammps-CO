"""
Master script: generates all 6 paper figures in sequence.
Run from any directory:
    python tools/plots/run_all_figures.py
Requirements: matplotlib, numpy, scipy
    pip install matplotlib numpy scipy
"""

import subprocess, sys, pathlib, time

ROOT = pathlib.Path(__file__).parent
scripts = [
    ROOT / "fig1_workflow.py",
    ROOT / "fig2_omp_coverage.py",
    ROOT / "fig3_build_flags.py",
    ROOT / "fig4_omp_scaling.py",
    ROOT / "fig5_a3_backport.py",
    ROOT / "fig6_summary.py",
]

def check_deps():
    missing = []
    for pkg in ["matplotlib", "numpy", "scipy"]:
        try:
            __import__(pkg)
        except ImportError:
            missing.append(pkg)
    if missing:
        print(f"Installing missing packages: {missing}")
        subprocess.check_call([sys.executable, "-m", "pip", "install", *missing])

check_deps()

total_start = time.time()
for i, script in enumerate(scripts, 1):
    print(f"\n[{i}/{len(scripts)}] Running {script.name} ...")
    t0 = time.time()
    result = subprocess.run([sys.executable, str(script)], capture_output=True, text=True)
    elapsed = time.time() - t0
    if result.returncode == 0:
        print(f"  OK  {result.stdout.strip()}  ({elapsed:.1f}s)")
    else:
        print(f"  FAIL  ERROR in {script.name}:")
        print(result.stderr[-600:])

total = time.time() - total_start
out_dir = ROOT
pdfs = list(out_dir.glob("figure*.pdf"))
pngs = list(out_dir.glob("figure*.png"))
print(f"\n{'='*55}")
print(f"All figures complete in {total:.1f}s")
print(f"  PDF files : {len(pdfs)}")
print(f"  PNG files : {len(pngs)}")
print(f"  Output dir: {out_dir}")
print(f"{'='*55}")
