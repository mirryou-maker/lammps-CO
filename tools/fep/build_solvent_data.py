"""
Build LAMMPS data files for FEP Li+ desolvation study.
Solvents: EC, PC, DME  (BTFE/FEC deferred - need F parameters)
Li+ ion: Joung & Cheatham (2008) SPC/E-optimized parameters

OPLS-AA parameters from:
  Jorgensen et al. JACS 1996 (core OPLS-AA)
  Canongia Lopes & Padua, JPC-B 2004 (carbonate/ether OPLS-AA)
  Masia et al. JPC-B 2004 (EC, PC specific)

Usage:
    python build_solvent_data.py --solvent ec  --nmol 200 --out data/ec_pure.data
    python build_solvent_data.py --solvent dme --nmol 200 --out data/dme_pure.data
"""

import numpy as np
import argparse, pathlib, textwrap, random, math

# ─────────────────────────────────────────────────────────────────────────────
# OPLS-AA parameter tables  (ε in kcal/mol, σ in Å, charge in e)
# ─────────────────────────────────────────────────────────────────────────────

# Atom-type definitions  key → (mass, eps, sig, charge, description)
ATOM_TYPES = {
    # Li+ Joung & Cheatham 2008 (SPC/E optimized, TIP3P compatible)
    "Li": (6.941,  0.3367,  1.4094, +1.00, "Li+ ion"),

    # EC / PC carbonyl carbon  (OPLS-AA #784 / Masia)
    "C_car": (12.011, 0.1050, 3.7500, +0.720, "Carbonyl C in cyclic carbonate"),
    # Carbonyl oxygen  (#466)
    "O_car": (15.999, 0.2100, 2.9600, -0.500, "Carbonyl O in cyclic carbonate"),
    # Ester (bridging) oxygen  (#467)
    "O_est": (15.999, 0.1700, 3.0000, -0.368, "Ester O in cyclic carbonate"),
    # Ring methylene CH2  (#145)
    "C_ec":  (12.011, 0.0660, 3.5000, +0.195, "Ring CH2 alpha to O in EC"),
    # Ring CH2 H  (#146)
    "H_ec":  ( 1.008, 0.0300, 2.5000, +0.000, "H on ring CH2 in EC"),

    # PC extra: methyl-bearing ring carbon
    "C_pc":  (12.011, 0.0660, 3.5000, +0.245, "Ring CH(CH3) in PC"),
    "C_me":  (12.011, 0.0660, 3.5000, -0.180, "Methyl C in PC"),
    "H_me":  ( 1.008, 0.0300, 2.5000, +0.060, "H on methyl in PC"),
    "H_pc":  ( 1.008, 0.0300, 2.5000, +0.000, "H on ring CH in PC"),

    # DME ether: CH3-O-CH2-CH2-O-CH3
    "CT":    (12.011, 0.0660, 3.5000, -0.100, "Methyl C in DME"),
    "H_CT":  ( 1.008, 0.0300, 2.5000, +0.030, "H on methyl in DME"),
    "OS":    (15.999, 0.1400, 3.0000, -0.400, "Ether O in DME"),
    "CT2":   (12.011, 0.0660, 3.5000, +0.250, "Methylene C alpha to O in DME"),
    "H_CT2": ( 1.008, 0.0300, 2.5000, +0.000, "H on CH2 in DME"),
}

# Bond parameters (r0 in Å, k in kcal/mol/Å²  OPLS-AA uses k/2*(r-r0)^2 → harmonic)
BOND_TYPES = {
    ("C_car","O_car"): (1.200, 1280.0),
    ("C_car","O_est"): (1.360,  900.0),
    ("O_est","C_ec"):  (1.430,  900.0),
    ("C_ec", "C_ec"):  (1.529,  900.0),
    ("C_ec", "H_ec"):  (1.090, 1404.0),
    # PC
    ("O_est","C_pc"):  (1.430,  900.0),
    ("C_pc", "C_me"):  (1.529,  900.0),
    ("C_pc", "H_pc"):  (1.090, 1404.0),
    ("C_me", "H_me"):  (1.090, 1404.0),
    ("C_pc", "C_ec"):  (1.529,  900.0),
    # DME
    ("CT",  "OS"):     (1.410,  900.0),
    ("CT",  "H_CT"):   (1.090, 1404.0),
    ("OS",  "CT2"):    (1.410,  900.0),
    ("CT2", "CT2"):    (1.529,  900.0),
    ("CT2", "H_CT2"):  (1.090, 1404.0),
}

# Angle parameters (theta0 in degrees, k in kcal/mol/rad²)
ANGLE_TYPES = {
    ("O_car","C_car","O_est"): (126.0, 160.0),
    ("C_car","O_est","C_ec"):  (109.0, 150.0),
    ("O_est","C_ec","C_ec"):   (109.5, 100.0),
    ("O_est","C_ec","H_ec"):   (109.5,  75.0),
    ("C_ec", "C_ec","H_ec"):   (109.5,  75.0),
    ("H_ec", "C_ec","H_ec"):   (107.8,  75.0),
    ("O_est","C_ec","O_est"):  (109.5, 100.0),  # symmetry
    # PC extras
    ("C_car","O_est","C_pc"):  (109.0, 150.0),
    ("O_est","C_pc","C_me"):   (109.5, 100.0),
    ("O_est","C_pc","H_pc"):   (109.5,  75.0),
    ("C_me", "C_pc","H_pc"):   (109.5,  75.0),
    ("C_pc", "C_me","H_me"):   (109.5,  75.0),
    ("H_me", "C_me","H_me"):   (107.8,  75.0),
    ("C_pc", "C_ec","H_ec"):   (109.5,  75.0),
    ("O_est","C_pc","C_ec"):   (109.5, 100.0),
    # DME
    ("CT", "OS","CT2"):        (109.5, 150.0),
    ("OS","CT","H_CT"):        (109.5,  75.0),
    ("OS","CT2","CT2"):        (109.5, 100.0),
    ("OS","CT2","H_CT2"):      (109.5,  75.0),
    ("CT2","CT2","H_CT2"):     (109.5,  75.0),
    ("H_CT","CT","H_CT"):      (107.8,  75.0),
    ("H_CT2","CT2","H_CT2"):   (107.8,  75.0),
    ("CT","CT2","H_CT2"):      (109.5,  75.0),  # unused but defensive
}

# ─────────────────────────────────────────────────────────────────────────────
# Molecular geometry builders  (returns list of (atype, x, y, z))
# ─────────────────────────────────────────────────────────────────────────────

def _rot(v, axis, theta):
    """Rodrigues rotation of vector v around axis by theta radians."""
    axis = axis / np.linalg.norm(axis)
    return (v * math.cos(theta)
            + np.cross(axis, v) * math.sin(theta)
            + axis * np.dot(axis, v) * (1 - math.cos(theta)))

def build_ec():
    """Ethylene carbonate: O=C1OCCO1 (5-membered ring + exo carbonyl O)"""
    # Ring: C_car - O_est - C_ec - C_ec - O_est (5-membered, planar)
    R = 1.50   # approximate ring radius for 5-membered ring
    theta0 = 2*math.pi / 5
    # Place ring atoms at vertices of regular pentagon
    ring_xyz = [(R*math.cos(i*theta0), R*math.sin(i*theta0), 0.0) for i in range(5)]
    # Atom order in ring: C_car(0), O_est(1), C_ec(2), C_ec(3), O_est(4)
    ring_types = ["C_car","O_est","C_ec","C_ec","O_est"]

    atoms = []
    for t, (x,y,z) in zip(ring_types, ring_xyz):
        atoms.append((t, x, y, z))

    # Exo carbonyl O: along z from C_car
    cx, cy, _ = ring_xyz[0]
    atoms.append(("O_car", cx, cy, 1.20))

    # H atoms on C_ec (atoms 2 and 3)
    for ci in [2, 3]:
        cx, cy, _ = ring_xyz[ci]
        # Two H above and below the ring plane
        atoms.append(("H_ec", cx + 0.50, cy, +0.89))
        atoms.append(("H_ec", cx - 0.50, cy, -0.89))

    return atoms   # 11 atoms total

def build_pc():
    """Propylene carbonate: O=C1OCC(C)O1 (5-membered ring, one CH2 → CHCH3)"""
    R = 1.50
    theta0 = 2*math.pi / 5
    ring_xyz = [(R*math.cos(i*theta0), R*math.sin(i*theta0), 0.0) for i in range(5)]
    ring_types = ["C_car","O_est","C_ec","C_pc","O_est"]

    atoms = []
    for t, (x,y,z) in zip(ring_types, ring_xyz):
        atoms.append((t, x, y, z))

    # Exo carbonyl O
    cx, cy, _ = ring_xyz[0]
    atoms.append(("O_car", cx, cy, 1.20))

    # C_ec (index 2): two H
    cx, cy, _ = ring_xyz[2]
    atoms.append(("H_ec", cx + 0.50, cy, +0.89))
    atoms.append(("H_ec", cx - 0.50, cy, -0.89))

    # C_pc (index 3): one H + methyl group
    cx, cy, _ = ring_xyz[3]
    atoms.append(("H_pc", cx, cy + 0.60, +0.89))
    # Methyl carbon (tetrahedral-ish placement)
    atoms.append(("C_me", cx + 1.54*math.cos(math.pi/6), cy - 1.54*math.sin(math.pi/6), 0.0))
    me_cx = cx + 1.54*math.cos(math.pi/6)
    me_cy = cy - 1.54*math.sin(math.pi/6)
    for ang in [0, 2*math.pi/3, 4*math.pi/3]:
        atoms.append(("H_me", me_cx + 1.09*math.cos(ang), me_cy + 1.09*math.sin(ang), 0.63))

    return atoms   # 15 atoms

def build_dme():
    """DME: CH3-O-CH2-CH2-O-CH3  (gauche conformation, C2h)"""
    # Place backbone: CT-OS-CT2-CT2-OS-CT along x-axis with standard bond angles
    # Bond lengths: C-O 1.41, C-C 1.53, C-H 1.09
    atoms = []
    # Backbone atoms (simplified all-trans)
    backbone = [
        ("CT",  -3.10, 0.0, 0.0),
        ("OS",  -1.69, 0.0, 0.0),
        ("CT2", -0.90, 1.11, 0.0),
        ("CT2",  0.54, 1.11, 0.0),
        ("OS",   1.33, 0.0, 0.0),
        ("CT",   2.74, 0.0, 0.0),
    ]
    atoms.extend(backbone)

    # H on first CT (methyl): 3 H in staggered positions
    for ang in [0, 2*math.pi/3, 4*math.pi/3]:
        atoms.append(("H_CT", -3.50, 1.03*math.cos(ang), 1.03*math.sin(ang)))

    # H on last CT (methyl)
    for ang in [0, 2*math.pi/3, 4*math.pi/3]:
        atoms.append(("H_CT", 3.14, 1.03*math.cos(ang), 1.03*math.sin(ang)))

    # H on CT2 (atom index 2, position (-0.90, 1.11, 0))
    atoms.append(("H_CT2", -1.30, 1.80, +0.89))
    atoms.append(("H_CT2", -1.30, 1.80, -0.89))

    # H on CT2 (atom index 3, position (0.54, 1.11, 0))
    atoms.append(("H_CT2", 0.94, 1.80, +0.89))
    atoms.append(("H_CT2", 0.94, 1.80, -0.89))

    return atoms   # 20 atoms

SOLVENT_BUILDERS = {
    "ec":  (build_ec,  "Ethylene carbonate"),
    "pc":  (build_pc,  "Propylene carbonate"),
    "dme": (build_dme, "1,2-Dimethoxyethane"),
}

# ─────────────────────────────────────────────────────────────────────────────
# Bond / angle topology builders  (returns list of (t1,t2,...) tuples for mol)
# ─────────────────────────────────────────────────────────────────────────────

def _bonds_angles_ec(offset):
    """Return (bonds, angles) for EC. offset = first atom index (1-based)."""
    o = offset
    # Atom mapping: 1=C_car, 2=O_est, 3=C_ec, 4=C_ec, 5=O_est, 6=O_car, 7-10=H_ec
    bonds = [
        (o+0, o+5, "C_car","O_est"),   # C_car-O_est(4)
        (o+0, o+1, "C_car","O_est"),   # C_car-O_est(1)
        (o+0, o+5, "C_car","O_car"),   # reuse key -> use index directly below
        (o+1, o+2, "O_est","C_ec"),
        (o+2, o+3, "C_ec","C_ec"),
        (o+3, o+4, "O_est","C_ec"),   # O_est(4)-C_ec(3) (reversed lookup ok)
        (o+2, o+6, "C_ec","H_ec"),
        (o+2, o+7, "C_ec","H_ec"),
        (o+3, o+8, "C_ec","H_ec"),
        (o+3, o+9, "C_ec","H_ec"),
    ]
    # Rebuild properly
    bonds = [
        (o+1, o+2, "C_car","O_est"),
        (o+1, o+5, "C_car","O_est"),
        (o+1, o+6, "C_car","O_car"),
        (o+2, o+3, "O_est","C_ec"),
        (o+3, o+4, "C_ec","C_ec"),
        (o+4, o+5, "O_est","C_ec"),
        (o+3, o+7, "C_ec","H_ec"),
        (o+3, o+8, "C_ec","H_ec"),
        (o+4, o+9, "C_ec","H_ec"),
        (o+4,o+10, "C_ec","H_ec"),
    ]
    angles = [
        (o+2, o+1, o+5, "O_est","C_car","O_est"),
        (o+2, o+1, o+6, "O_est","C_car","O_car"),
        (o+5, o+1, o+6, "O_est","C_car","O_car"),
        (o+1, o+2, o+3, "C_car","O_est","C_ec"),
        (o+2, o+3, o+4, "O_est","C_ec","C_ec"),
        (o+2, o+3, o+7, "O_est","C_ec","H_ec"),
        (o+2, o+3, o+8, "O_est","C_ec","H_ec"),
        (o+4, o+3, o+7, "C_ec","C_ec","H_ec"),
        (o+4, o+3, o+8, "C_ec","C_ec","H_ec"),
        (o+7, o+3, o+8, "H_ec","C_ec","H_ec"),
        (o+3, o+4, o+5, "C_ec","C_ec","O_est"),
        (o+3, o+4, o+9, "C_ec","C_ec","H_ec"),
        (o+3, o+4,o+10, "C_ec","C_ec","H_ec"),
        (o+5, o+4, o+9, "O_est","C_ec","H_ec"),
        (o+5, o+4,o+10, "O_est","C_ec","H_ec"),
        (o+9, o+4,o+10, "H_ec","C_ec","H_ec"),
        (o+1, o+5, o+4, "C_car","O_est","C_ec"),
    ]
    return bonds, angles

# ─────────────────────────────────────────────────────────────────────────────
# Simple random packing (avoid close contacts)
# ─────────────────────────────────────────────────────────────────────────────

def random_pack(mol_atoms, n_mol, box, n_li=1, seed=42, min_dist=3.5):
    """
    Randomly place n_mol copies of mol_atoms + n_li Li+ ions in box.
    box = (lx, ly, lz)
    Returns list of (atype, x, y, z) for all atoms.
    """
    rng = random.Random(seed)
    lx, ly, lz = box
    mol_arr = np.array([(x,y,z) for _,x,y,z in mol_atoms])
    mol_types = [t for t,_,_,_ in mol_atoms]

    placed_com = []   # list of (cx,cy,cz) for placed molecules
    all_atoms = []

    max_attempts = 5000

    for m in range(n_mol):
        for attempt in range(max_attempts):
            # Random translation
            tx = rng.uniform(0, lx)
            ty = rng.uniform(0, ly)
            tz = rng.uniform(0, lz)
            # Random rotation (random quaternion)
            u1, u2, u3 = rng.random(), rng.random(), rng.random()
            q = np.array([
                math.sqrt(1-u1)*math.sin(2*math.pi*u2),
                math.sqrt(1-u1)*math.cos(2*math.pi*u2),
                math.sqrt(u1)*math.sin(2*math.pi*u3),
                math.sqrt(u1)*math.cos(2*math.pi*u3),
            ])
            # Rotation matrix from quaternion
            w,x,y,z2 = q
            R = np.array([
                [1-2*(y*y+z2*z2), 2*(x*y-w*z2),   2*(x*z2+w*y)],
                [2*(x*y+w*z2),   1-2*(x*x+z2*z2), 2*(y*z2-w*x)],
                [2*(x*z2-w*y),   2*(y*z2+w*x),   1-2*(x*x+y*y)],
            ])
            rotated = mol_arr @ R.T + np.array([tx, ty, tz])

            # Check min distance to already placed molecule COMs
            ok = True
            com = rotated.mean(axis=0)
            for pc in placed_com:
                d = np.linalg.norm(com - np.array(pc))
                if d < min_dist * 1.5:
                    ok = False
                    break
            if ok:
                placed_com.append(tuple(com))
                for t, (ax,ay,az) in zip(mol_types, rotated):
                    all_atoms.append((t, ax % lx, ay % ly, az % lz))
                break
        else:
            # Give up trying to avoid overlaps — just place randomly
            tx = rng.uniform(0, lx)
            ty = rng.uniform(0, ly)
            tz = rng.uniform(0, lz)
            for t, (ax,ay,az) in zip(mol_types, mol_arr):
                all_atoms.append((t, (ax+tx)%lx, (ay+ty)%ly, (az+tz)%lz))

    # Add Li+ ions
    for _ in range(n_li):
        all_atoms.append(("Li", rng.uniform(0,lx), rng.uniform(0,ly), rng.uniform(0,lz)))

    return all_atoms

# ─────────────────────────────────────────────────────────────────────────────
# LAMMPS data file writer
# ─────────────────────────────────────────────────────────────────────────────

def write_lammps_data(filename, all_atoms, mol_atoms, n_mol, box, solvent_name):
    """Write a full LAMMPS data file with atom types, masses, pair coeffs."""
    lx, ly, lz = box
    n_atoms_total = len(all_atoms)

    # Build type index
    type_names_in_mol = list(dict.fromkeys(t for t,_,_,_ in mol_atoms))
    if "Li" not in type_names_in_mol:
        type_names_in_mol.append("Li")
    type_idx = {t: i+1 for i, t in enumerate(type_names_in_mol)}
    n_types = len(type_names_in_mol)

    # Atom count per molecule for molecule ID assignment
    natoms_per_mol = len(mol_atoms)

    with open(filename, "w") as f:
        f.write(f"# LAMMPS data file: {solvent_name} + Li+ ({n_mol} mol + 1 Li+)\n")
        f.write(f"# Generated by build_solvent_data.py\n\n")
        f.write(f"{n_atoms_total} atoms\n")
        f.write(f"{n_types} atom types\n\n")
        f.write(f"0.0 {lx:.4f} xlo xhi\n")
        f.write(f"0.0 {ly:.4f} ylo yhi\n")
        f.write(f"0.0 {lz:.4f} zlo zhi\n\n")

        f.write("Masses\n\n")
        for t, idx in type_idx.items():
            mass = ATOM_TYPES[t][0]
            f.write(f"  {idx}  {mass:.4f}  # {t}\n")
        f.write("\n")

        f.write("Pair Coeffs  # lj/cut/soft  (eps kcal/mol, sigma Ang, lambda)\n\n")
        for t, idx in type_idx.items():
            _, eps, sig, _, desc = ATOM_TYPES[t]
            f.write(f"  {idx}  {eps:.6f}  {sig:.4f}  # {t}: {desc}\n")
        f.write("\n")

        f.write("Atoms  # full\n\n")
        for i, (atype, x, y, z) in enumerate(all_atoms, 1):
            mol_id = ((i-1) // natoms_per_mol) + 1 if atype != "Li" else n_mol + 1
            if atype == "Li":
                mol_id = n_mol + 1
            tidx = type_idx[atype]
            charge = ATOM_TYPES[atype][3]
            f.write(f"  {i}  {mol_id}  {tidx}  {charge:.4f}  {x:.6f}  {y:.6f}  {z:.6f}  # {atype}\n")

    print(f"Written: {filename}  ({n_atoms_total} atoms, {n_types} types)")
    return type_idx, type_names_in_mol

# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--solvent", choices=["ec","pc","dme"], required=True)
    parser.add_argument("--nmol", type=int, default=200)
    parser.add_argument("--out", required=True)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    builder, desc = SOLVENT_BUILDERS[args.solvent]
    mol_atoms = builder()
    print(f"Molecule: {desc}  ({len(mol_atoms)} atoms/molecule)")

    # Box size: target density ~1.1 g/cm3 for carbonates, ~0.87 for DME
    # Volume per molecule = M / (rho * NA)
    mol_mass = sum(ATOM_TYPES[t][0] for t,_,_,_ in mol_atoms)
    rho = 0.87e-24 if args.solvent == "dme" else 1.10e-24  # g/atom
    rho_mol = rho * 6.022e23  # g/mol per atom → convert to g/cm3 per molecule
    vol_mol = mol_mass / (rho_mol * 6.022e23) * 1e24  # cm3 → Å3
    box_vol = vol_mol * args.nmol
    box_len = box_vol ** (1/3)
    box = (box_len, box_len, box_len)
    print(f"Box: {box_len:.2f} Å  (target ρ for {args.solvent.upper()}, {args.nmol} mol)")

    print(f"Packing {args.nmol} molecules + 1 Li+  (seed={args.seed}) ...")
    all_atoms = random_pack(mol_atoms, args.nmol, box, n_li=1, seed=args.seed)

    pathlib.Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    write_lammps_data(args.out, all_atoms, mol_atoms, args.nmol, box, desc)

if __name__ == "__main__":
    main()
