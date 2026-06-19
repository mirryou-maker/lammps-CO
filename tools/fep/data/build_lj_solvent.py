"""
Build a simple LJ solvent + Li+ system for FEP pipeline validation.
Each "solvent bead" represents one EC-like molecule (united-atom CG).
No bonds/angles needed — pure pair interactions.

EC CG bead: sigma=5.0 Ang, eps=0.5 kcal/mol, charge=0 (non-polar CG)
Li+ ion:    sigma=1.4094 Ang, eps=0.3367 kcal/mol, charge=+1

The FEP measures ΔG for decoupling Li+ from CG solvent.
Purpose: validate the 20-window FEP OMP speedup pipeline before running
         full all-atom OPLS-AA simulations.
"""
import numpy as np, pathlib, math, random

def build_lj_fep_data(n_solvent=500, out="lj_solvent_li.data", seed=42, rho=0.87):
    """
    n_solvent: number of CG solvent beads
    rho: number density in Ang^-3 (0.87/M_EC*NA ≈ 0.0040 Ang^-3 for EC)
    """
    rng = random.Random(seed)
    # Target number density for liquid EC: ~1.32 g/cm3, M=88 g/mol
    # n/V = rho_mass * NA / M = 1.32e-24 * 6.022e23 / 88 = 9.03e-3 per Ang3
    # But for CG sigma=5.0, reduced density rho*sigma^3 ~ 0.5 → liquid
    sigma_solv = 5.0  # Ang
    rho_reduced = 0.50  # liquid-like LJ reduced density
    n_density = rho_reduced / (sigma_solv**3)  # per Ang^3
    vol = n_solvent / n_density
    L = vol**(1/3)

    # Random placement with simple rejection (min_dist = sigma_solv * 0.8)
    min_d = sigma_solv * 0.85
    positions = []
    for _ in range(n_solvent):
        for attempt in range(10000):
            pos = np.array([rng.uniform(0,L), rng.uniform(0,L), rng.uniform(0,L)])
            ok = True
            for prev in positions:
                d = np.linalg.norm(pos - prev)
                if d < min_d: ok = False; break
            if ok: positions.append(pos); break
        else:
            positions.append(np.array([rng.uniform(0,L), rng.uniform(0,L), rng.uniform(0,L)]))

    # Li+ at box center
    li_pos = np.array([L/2, L/2, L/2])

    total = n_solvent + 1
    with open(out, "w") as f:
        f.write(f"# LJ CG solvent + Li+  ({n_solvent} beads + 1 Li+)\n")
        f.write(f"# Box length: {L:.3f} Ang, rho_red={rho_reduced}\n\n")
        f.write(f"{total} atoms\n")
        f.write(f"2 atom types\n\n")
        f.write(f"0.0 {L:.4f} xlo xhi\n")
        f.write(f"0.0 {L:.4f} ylo yhi\n")
        f.write(f"0.0 {L:.4f} zlo zhi\n\n")
        f.write("Masses\n\n")
        f.write(f"  1  88.062  # EC_CG (ethylene carbonate CG bead)\n")
        f.write(f"  2   6.941  # Li+\n\n")
        f.write("Atoms  # charge\n\n")
        for i, pos in enumerate(positions, 1):
            f.write(f"  {i}  1  0.000  {pos[0]:.6f}  {pos[1]:.6f}  {pos[2]:.6f}\n")
        # Li+
        f.write(f"  {total}  2  1.000  {li_pos[0]:.6f}  {li_pos[1]:.6f}  {li_pos[2]:.6f}\n")

    print(f"Written {out}: {n_solvent} CG beads + 1 Li+, box={L:.2f} Ang")
    return L

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=500)
    p.add_argument("--out", default="lj_solvent_li.data")
    args = p.parse_args()
    build_lj_fep_data(n_solvent=args.n, out=args.out)
