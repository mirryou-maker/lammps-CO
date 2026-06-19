"""
build_opls_systems.py  — OPLS-AA all-atom LAMMPS data files for FEP.

Generates: ec_li_aa.data, pc_li_aa.data, dme_li_aa.data
Each contains N solvent molecules + 1 Li+ ion.

Force field: OPLS-AA (Jorgensen 1996) + Joung & Cheatham (2008) for Li+.
Pair style:  lj/cut/coul/long/soft 2 0.5 0.3 10.0  (FEP ready)
"""

import numpy as np
import pathlib
import math

RNG = np.random.default_rng(42)

# ─── OPLS-AA single-site LJ parameters ─────────────────────────────────────
# (mass, sigma_Ang, eps_kcalmol)
LJ = {
    'Li':    (6.941,  1.4094, 0.33670),   # Joung & Cheatham 2008
    'C_co3': (12.011, 3.7500, 0.10500),   # carbonyl C in cyclic carbonate
    'O_co3': (15.999, 2.9600, 0.21000),   # exocyclic C=O oxygen
    'O_ete': (15.999, 3.0000, 0.17000),   # ring/ether O
    'C_sp3': (12.011, 3.5000, 0.06600),   # generic sp3 C (CH2, CH3, CH)
    'H_c':   (1.008,  2.5000, 0.03000),   # H on C
}

# ─── Bond parameters  K (kcal/mol/Ang^2), r0 (Ang) ─────────────────────────
BOND_K = {
    ('C_co3','O_co3'): (570.0, 1.218),
    ('C_co3','O_ete'): (450.0, 1.370),
    ('C_sp3','O_ete'): (320.0, 1.410),
    ('C_sp3','C_sp3'): (224.0, 1.529),
    ('C_sp3','H_c'):   (340.0, 1.090),
}

# ─── Angle parameters  K (kcal/mol/rad^2), theta0 (deg) ────────────────────
ANGL_K = {
    ('O_co3','C_co3','O_ete'): (85.0, 125.0),
    ('O_ete','C_co3','O_ete'): (85.0, 109.5),
    ('C_co3','O_ete','C_sp3'): (60.0, 109.5),
    ('O_ete','C_sp3','C_sp3'): (50.0, 107.0),
    ('O_ete','C_sp3','H_c'):   (50.0, 110.0),
    ('C_sp3','C_sp3','O_ete'): (50.0, 109.5),
    ('C_sp3','C_sp3','H_c'):   (37.5, 110.7),
    ('H_c',  'C_sp3','H_c'):   (33.0, 107.8),
    ('C_sp3','O_ete','C_sp3'): (60.0, 111.8),   # DME ether angle
    ('H_c',  'C_sp3','C_sp3'): (37.5, 110.7),   # alias
    ('H_c',  'C_sp3','O_ete'): (50.0, 110.0),   # alias
}

# ─── Dihedral parameters  OPLS style: E = c1+c2+c3+c4 ──────────────────────
# Key = sorted(center two atoms); c1..c4 in kcal/mol
DIHE_K = {
    ('C_co3','O_ete'): (0.0,  0.0,  0.0,  0.0),   # ring constraint
    ('O_ete','C_sp3'): (0.0,  0.0,  0.35, 0.0),   # H-C-O-C ether
    ('C_sp3','C_sp3'): (0.0,  0.0,  0.30, 0.0),   # C-C torsion
    ('C_co3','O_co3'): (0.0,  2.70, 0.0,  0.0),   # exo O planarity
    ('C_sp3','O_co3'): (0.0,  0.0,  0.0,  0.0),   # not used
}

# ─── Geometry helpers ────────────────────────────────────────────────────────

def _add_methyl_H(c_pos, bond_vec):
    """Return 3 H positions in tetrahedral arrangement from methyl C.
    bond_vec: unit vector from CH3-C toward its bonded neighbor.
    """
    bond_vec = np.asarray(bond_vec, float)
    bond_vec /= np.linalg.norm(bond_vec)
    perp = np.cross(bond_vec, [0, 1, 0])
    if np.linalg.norm(perp) < 0.1:
        perp = np.cross(bond_vec, [1, 0, 0])
    perp /= np.linalg.norm(perp)
    perp2 = np.cross(bond_vec, perp)
    theta = math.radians(109.5)
    H = []
    for phi_deg in [0, 120, 240]:
        phi = math.radians(phi_deg)
        h_dir = (math.cos(theta) * (-bond_vec) +
                 math.sin(theta) * (math.cos(phi) * perp + math.sin(phi) * perp2))
        H.append(c_pos + 1.09 * h_dir)
    return H


def _rot_rand(pos_list):
    """Apply a random rotation to a list of positions (Nx3 array)."""
    u = RNG.standard_normal(3); u /= np.linalg.norm(u)
    v = RNG.standard_normal(3); v -= np.dot(v, u) * u; v /= np.linalg.norm(v)
    w = np.cross(u, v)
    R = np.column_stack([u, v, w])
    return [R @ p for p in pos_list]


# ─── Molecule templates ──────────────────────────────────────────────────────

def ec_template():
    """Ethylene carbonate (EC, C3H4O3) OPLS-AA.

    Charges sum to 0.00:
      C_co3:+0.76, O_co3:-0.50, O_ete×2:-0.35, C_sp3×2:+0.10, H_c×4:+0.06
    """
    R = 1.216   # regular-pentagon circumradius (side ~1.43 Ang)
    def p(deg): return R * np.array([math.cos(math.radians(deg)),
                                     math.sin(math.radians(deg)), 0.0])
    C1  = p(90)                              # carbonyl C
    O2  = p(162)                             # ring O (left)
    C4  = p(234)                             # ring CH2 (left)
    C5  = p(306)                             # ring CH2 (right)
    O3  = p(18)                              # ring O (right)
    O1  = C1 + 1.22 * C1 / np.linalg.norm(C1)  # exo C=O

    H41 = C4 + np.array([0., 0., +1.09])
    H42 = C4 + np.array([0., 0., -1.09])
    H51 = C5 + np.array([0., 0., +1.09])
    H52 = C5 + np.array([0., 0., -1.09])

    atoms = [
        # (type, charge, xyz)          idx
        ('C_co3', +0.76, C1),   # 0
        ('O_co3', -0.50, O1),   # 1
        ('O_ete', -0.35, O2),   # 2
        ('O_ete', -0.35, O3),   # 3
        ('C_sp3', +0.10, C4),   # 4
        ('C_sp3', +0.10, C5),   # 5
        ('H_c',   +0.06, H41),  # 6
        ('H_c',   +0.06, H42),  # 7
        ('H_c',   +0.06, H51),  # 8
        ('H_c',   +0.06, H52),  # 9
    ]
    assert abs(sum(q for _,q,_ in atoms)) < 1e-9, "EC charge not neutral"

    # bonds: (i, j) 0-based within molecule
    bonds = [(0,1),(0,2),(0,3),(2,4),(3,5),(4,5),(4,6),(4,7),(5,8),(5,9)]

    # angles: (i, j_center, k)
    angles = [
        (1,0,2),(1,0,3),(2,0,3),
        (0,2,4),(0,3,5),
        (2,4,5),(3,5,4),
        (2,4,6),(2,4,7),(3,5,8),(3,5,9),
        (5,4,6),(5,4,7),(4,5,8),(4,5,9),
        (6,4,7),(8,5,9),
    ]

    # dihedrals: (i,j,k,l)
    dihedrals = [
        (1,0,2,4),(1,0,3,5),           # O_exo–C–O_ring–C (planarity)
        (3,0,2,4),(2,0,3,5),           # O_ring–C–O_ring–C (ring)
        (0,2,4,5),(0,3,5,4),           # C–O_ring–CH2–CH2
        (0,2,4,6),(0,2,4,7),           # C–O–CH2–H
        (0,3,5,8),(0,3,5,9),
        (2,4,5,3),(2,4,5,8),(2,4,5,9), # C–CH2–CH2–O and C–CH2–CH2–H
        (3,5,4,2),(3,5,4,6),(3,5,4,7),
        (6,4,5,8),(6,4,5,9),(7,4,5,8),(7,4,5,9),
    ]
    return atoms, bonds, angles, dihedrals


def pc_template():
    """Propylene carbonate (PC, C4H6O3) OPLS-AA.

    EC ring but C5 → C_CHr (1 H) with extra CH3 substituent.
    Charges: same ring as EC; CHr:+0.02, H_CHr:+0.06, C_CH3:+0.02, H_CH3×3:+0.04
    Total: 0.76-0.50-0.35-0.35+0.10+0.02+0.06+0.06+0.06+0.02+0.04×3 = 0.00 ✓
    """
    R = 1.216
    def p(deg): return R * np.array([math.cos(math.radians(deg)),
                                     math.sin(math.radians(deg)), 0.0])
    C1  = p(90)
    O2  = p(162)
    C4  = p(234)                        # ring CH2 (left)
    C5  = p(306)                        # ring CH (right, bearing CH3)
    O3  = p(18)
    O1  = C1 + 1.22 * C1 / np.linalg.norm(C1)

    H41 = C4 + np.array([0., 0., +1.09])
    H42 = C4 + np.array([0., 0., -1.09])
    H51 = C5 + np.array([0., 0., +1.09])   # single H on CHr

    # CH3 substituent on C5: roughly in the -z direction from C5
    C_me = C5 + 1.54 * np.array([0.0, -0.374, -0.928])
    C_me /= 1.0   # keep position
    # Normalize direction for H placement
    bond_to_C5 = (C5 - C_me); bond_to_C5 /= np.linalg.norm(bond_to_C5)
    H_me = _add_methyl_H(C_me, bond_to_C5)

    atoms = [
        ('C_co3', +0.76, C1),    # 0
        ('O_co3', -0.50, O1),    # 1
        ('O_ete', -0.35, O2),    # 2
        ('O_ete', -0.35, O3),    # 3
        ('C_sp3', +0.10, C4),    # 4   ring CH2
        ('C_sp3', +0.02, C5),    # 5   ring CHr (bearing CH3)
        ('H_c',   +0.06, H41),   # 6
        ('H_c',   +0.06, H42),   # 7
        ('H_c',   +0.06, H51),   # 8   H on CHr
        ('C_sp3', +0.02, C_me),  # 9   methyl C
        ('H_c',   +0.04, H_me[0]),# 10
        ('H_c',   +0.04, H_me[1]),# 11
        ('H_c',   +0.04, H_me[2]),# 12
    ]
    assert abs(sum(q for _,q,_ in atoms)) < 1e-9, f"PC charge: {sum(q for _,q,_ in atoms)}"

    bonds = [
        (0,1),(0,2),(0,3),(2,4),(3,5),(4,5),
        (4,6),(4,7),(5,8),(5,9),
        (9,10),(9,11),(9,12),
    ]
    angles = [
        (1,0,2),(1,0,3),(2,0,3),
        (0,2,4),(0,3,5),
        (2,4,5),(3,5,4),
        (2,4,6),(2,4,7),(3,5,8),(3,5,9),
        (5,4,6),(5,4,7),(4,5,8),(4,5,9),
        (6,4,7),
        (5,9,10),(5,9,11),(5,9,12),
        (10,9,11),(10,9,12),(11,9,12),
    ]
    dihedrals = [
        (1,0,2,4),(1,0,3,5),
        (3,0,2,4),(2,0,3,5),
        (0,2,4,5),(0,3,5,4),
        (0,2,4,6),(0,2,4,7),(0,3,5,8),(0,3,5,9),
        (2,4,5,3),(2,4,5,8),(2,4,5,9),
        (3,5,4,2),(3,5,4,6),(3,5,4,7),
        (6,4,5,8),(6,4,5,9),(7,4,5,8),(7,4,5,9),
        (4,5,9,10),(4,5,9,11),(4,5,9,12),
        (8,5,9,10),(8,5,9,11),(8,5,9,12),
        (3,5,9,10),(3,5,9,11),(3,5,9,12),
    ]
    return atoms, bonds, angles, dihedrals


def dme_template():
    """1,2-dimethoxyethane (DME, C4H10O2) OPLS-AA, all-trans.

    Charges:
      C_sp3(CH3)×2:+0.022, O_ete×2:-0.400, C_sp3(CH2)×2:+0.164,
      H_c(CH3)×6:+0.034, H_c(CH2)×4:+0.056
    Total = 0.044-0.800+0.328+0.204+0.224 = 0.000 ✓
    """
    CT_a = np.array([-2.820, 0., 0.])
    O_a  = np.array([-1.410, 0., 0.])
    CC_a = np.array([ 0.000, 0., 0.])
    CC_b = np.array([ 1.529, 0., 0.])
    O_b  = np.array([ 2.939, 0., 0.])
    CT_b = np.array([ 4.349, 0., 0.])

    # H on left CH3: bond direction CT_a→O_a = +x
    H_CT_a = _add_methyl_H(CT_a, np.array([+1., 0., 0.]))
    # H on left CH2 (2 H, tetrahedral around x-bond)
    H_CC_a = [CC_a + np.array([0., +0.890, +0.580]),
              CC_a + np.array([0., -0.890, +0.580])]
    # H on right CH2
    H_CC_b = [CC_b + np.array([0., +0.890, -0.580]),
              CC_b + np.array([0., -0.890, -0.580])]
    # H on right CH3: bond direction CT_b→O_b = -x
    H_CT_b = _add_methyl_H(CT_b, np.array([-1., 0., 0.]))

    atoms = [
        ('C_sp3', +0.022, CT_a),    # 0  left CH3
        ('O_ete', -0.400, O_a),     # 1  left O
        ('C_sp3', +0.164, CC_a),    # 2  left CH2
        ('C_sp3', +0.164, CC_b),    # 3  right CH2
        ('O_ete', -0.400, O_b),     # 4  right O
        ('C_sp3', +0.022, CT_b),    # 5  right CH3
        ('H_c',   +0.034, H_CT_a[0]),# 6
        ('H_c',   +0.034, H_CT_a[1]),# 7
        ('H_c',   +0.034, H_CT_a[2]),# 8
        ('H_c',   +0.056, H_CC_a[0]),# 9
        ('H_c',   +0.056, H_CC_a[1]),# 10
        ('H_c',   +0.056, H_CC_b[0]),# 11
        ('H_c',   +0.056, H_CC_b[1]),# 12
        ('H_c',   +0.034, H_CT_b[0]),# 13
        ('H_c',   +0.034, H_CT_b[1]),# 14
        ('H_c',   +0.034, H_CT_b[2]),# 15
    ]
    assert abs(sum(q for _,q,_ in atoms)) < 1e-9, "DME charge not neutral"

    bonds = [
        (0,1),(1,2),(2,3),(3,4),(4,5),
        (0,6),(0,7),(0,8),
        (2,9),(2,10),
        (3,11),(3,12),
        (5,13),(5,14),(5,15),
    ]
    angles = [
        (0,1,2),(1,2,3),(2,3,4),(3,4,5),      # backbone
        (1,0,6),(1,0,7),(1,0,8),               # O-CH3-H
        (6,0,7),(6,0,8),(7,0,8),               # H-CH3-H
        (0,1,2),(1,2,3),                        # C-O-C
        (1,2,9),(1,2,10),(1,2,3),              # O-CH2-H, O-CH2-C
        (3,2,9),(3,2,10),                      # C-CH2-H (right side)
        (9,2,10),                              # H-CH2-H
        (4,3,2),(4,3,11),(4,3,12),(2,3,11),(2,3,12),(11,3,12),
        (3,4,5),(4,5,13),(4,5,14),(4,5,15),
        (13,5,14),(13,5,15),(14,5,15),
    ]
    dihedrals = [
        (0,1,2,3),(0,1,2,9),(0,1,2,10),        # C-O-CH2-C, C-O-CH2-H
        (1,2,3,4),(1,2,3,11),(1,2,3,12),       # O-CH2-CH2-O, O-CH2-CH2-H
        (2,3,4,5),(2,3,4,11),(2,3,4,12),       # CH2-CH2-O-CH3, CH2-CH2-O-H
        (9,2,3,4),(9,2,3,11),(9,2,3,12),       # H-CH2-CH2-O
        (10,2,3,4),(10,2,3,11),(10,2,3,12),
        (3,4,5,13),(3,4,5,14),(3,4,5,15),      # C-O-CH3-H
        (6,0,1,2),(7,0,1,2),(8,0,1,2),         # H-CH3-O-C
        (9,2,1,0),(10,2,1,0),                  # H-CH2-O-C
    ]
    return atoms, bonds, angles, dihedrals


# ─── Topology → index maps ───────────────────────────────────────────────────

def _bond_type_key(at1, at2):
    return tuple(sorted([at1, at2]))

def _angle_type_key(at1, at2, at3):
    if at1 > at3: at1, at3 = at3, at1
    return (at1, at2, at3)

def _dihedral_type_key(at1, at2, at3, at4):
    return tuple(sorted([at2, at3]))


# ─── System builder ─────────────────────────────────────────────────────────

def build_system(mol_fn, n_mol, out_path, solvent_name, L_scale=1.5):
    """Pack n_mol molecules + 1 Li+ into a box; write LAMMPS data."""

    atoms0, bonds0, angles0, dihedrals0 = mol_fn()
    n_at = len(atoms0)

    # Collect unique FF types
    atom_types_all = [LJ['Li']] + [LJ[at] for at, _, _ in atoms0]
    atype_names = ['Li'] + [at for at, _, _ in atoms0]
    unique_atypes = list(dict.fromkeys(atype_names[1:]))  # preserve order, no Li yet
    atypes_ordered = ['Li'] + unique_atypes

    atype_idx = {name: i+1 for i, name in enumerate(atypes_ordered)}

    # Bond/angle/dihedral types from one molecule
    def bond_params(i, j):
        k = _bond_type_key(atoms0[i][0], atoms0[j][0])
        return BOND_K.get(k, BOND_K.get((k[1],k[0]), (200.0, 1.40)))

    def angle_params(i, j, k):
        key = _angle_type_key(atoms0[i][0], atoms0[j][0], atoms0[k][0])
        v = ANGL_K.get(key)
        if v is None:
            key2 = (key[2], key[1], key[0])
            v = ANGL_K.get(key2, (50.0, 109.5))
        return v

    def dihe_params(i, j, k, l):
        key = _dihedral_type_key(atoms0[i][0], atoms0[j][0], atoms0[k][0], atoms0[l][0])
        return DIHE_K.get(key, (0.0, 0.0, 0.3, 0.0))

    # Enumerate unique bond/angle/dihedral types
    btypes = {}
    for (i,j) in bonds0:
        k = _bond_type_key(atoms0[i][0], atoms0[j][0])
        if k not in btypes:
            btypes[k] = (len(btypes)+1, bond_params(i,j))

    atypes_ang = {}
    for (i,j,k) in angles0:
        key = _angle_type_key(atoms0[i][0], atoms0[j][0], atoms0[k][0])
        if key not in atypes_ang:
            atypes_ang[key] = (len(atypes_ang)+1, angle_params(i,j,k))

    dtypes = {}
    for (i,j,k,l) in dihedrals0:
        key = _dihedral_type_key(atoms0[i][0], atoms0[j][0], atoms0[k][0], atoms0[l][0])
        if key not in dtypes:
            dtypes[key] = (len(dtypes)+1, dihe_params(i,j,k,l))

    # ── Compute box size ──────────────────────────────────────────────────
    # Target density based on solvent
    densities = {'EC': 1.32, 'DME': 0.867, 'PC': 1.195}
    m_mol = {'EC': 88.062, 'DME': 90.122, 'PC': 102.089}
    rho = densities[solvent_name]
    M   = m_mol[solvent_name]
    NA  = 6.02214e23
    V_target = (n_mol * M / (NA * rho * 1e-24))  # Ang^3
    L_target  = V_target ** (1/3)
    L_box     = L_target * L_scale   # inflate for easier packing

    # ── Place molecules on grid ───────────────────────────────────────────
    grid_n = math.ceil(n_mol ** (1/3)) + 1
    spacing = L_box / grid_n
    grid_pts = [(ix*spacing, iy*spacing, iz*spacing)
                for ix in range(grid_n)
                for iy in range(grid_n)
                for iz in range(grid_n)][:n_mol]

    # Center molecule at origin
    com = np.mean([pos for _, _, pos in atoms0], axis=0)
    mol_centered = [(at, q, pos - com) for at, q, pos in atoms0]

    all_atoms  = []  # (mol_id, atype_idx, charge, x, y, z)
    all_bonds  = []  # (btype_idx, i_global, j_global)
    all_angles = []  # (atype_idx, i, j, k)
    all_dihe   = []  # (dtype_idx, i, j, k, l)

    for mol_id, (gx, gy, gz) in enumerate(grid_pts, start=1):
        center = np.array([gx, gy, gz])
        rotated = _rot_rand([pos for _, _, pos in mol_centered])
        # small random jitter
        jitter = RNG.uniform(-0.5, 0.5, 3)
        base_idx = (mol_id - 1) * n_at + 1  # 1-based global atom index

        for local_i, (at, q, _) in enumerate(mol_centered):
            pos = rotated[local_i] + center + jitter
            # wrap into box
            pos = pos % L_box
            all_atoms.append((mol_id, atype_idx[at], q, *pos))

        for (i, j) in bonds0:
            k = _bond_type_key(mol_centered[i][0], mol_centered[j][0])
            btype = btypes[k][0]
            all_bonds.append((btype, base_idx+i, base_idx+j))

        for (i, j, k) in angles0:
            key = _angle_type_key(mol_centered[i][0], mol_centered[j][0], mol_centered[k][0])
            v = atypes_ang.get(key)
            if v is None:
                key2 = (key[2], key[1], key[0])
                v = atypes_ang.get(key2, (1, (50.0,109.5)))
            atype_a = v[0]
            all_angles.append((atype_a, base_idx+i, base_idx+j, base_idx+k))

        for (i, j, k, l) in dihedrals0:
            key = _dihedral_type_key(mol_centered[i][0], mol_centered[j][0],
                                      mol_centered[k][0], mol_centered[l][0])
            dtype = dtypes[key][0]
            all_dihe.append((dtype, base_idx+i, base_idx+j, base_idx+k, base_idx+l))

    # Li+ at box center
    li_mol_id = n_mol + 1
    li_pos    = (L_box/2, L_box/2, L_box/2)
    all_atoms.append((li_mol_id, atype_idx['Li'], +1.0, *li_pos))

    n_atoms_total = len(all_atoms)
    n_bonds_total = len(all_bonds)
    n_angl_total  = len(all_angles)
    n_dihe_total  = len(all_dihe)

    # ── Compute pair coefficients (LB mixing for Li+ cross-pairs) ─────────
    # All pairs (Li type 1 vs solvent types)
    def lb_mix(name_a, name_b):
        _, sa, ea = LJ[name_a]
        _, sb, eb = LJ[name_b]
        return math.sqrt(ea*eb), (sa+sb)/2.0

    # ── Write LAMMPS data file ────────────────────────────────────────────
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        f.write(f"# OPLS-AA {solvent_name}+Li+ system  "
                f"({n_mol} molecules + 1 Li+, L={L_box:.2f} Ang)\n")
        f.write(f"# FEP: annihilate Li+-solvent coupling with lj/cut/coul/long/soft\n\n")
        f.write(f"{n_atoms_total} atoms\n")
        f.write(f"{n_bonds_total} bonds\n")
        f.write(f"{n_angl_total} angles\n")
        f.write(f"{n_dihe_total} dihedrals\n")
        f.write(f"0 impropers\n\n")
        f.write(f"{len(atypes_ordered)} atom types\n")
        f.write(f"{len(btypes)} bond types\n")
        f.write(f"{len(atypes_ang)} angle types\n")
        f.write(f"{len(dtypes)} dihedral types\n\n")
        f.write(f"0.000 {L_box:.4f} xlo xhi\n")
        f.write(f"0.000 {L_box:.4f} ylo yhi\n")
        f.write(f"0.000 {L_box:.4f} zlo zhi\n\n")

        f.write("Masses\n\n")
        for name in atypes_ordered:
            mass, _, _ = LJ[name]
            f.write(f"  {atype_idx[name]}  {mass:.3f}  # {name}\n")
        f.write("\n")

        f.write("Pair Coeffs  # lj/cut/coul/long/soft: epsilon sigma lambda\n\n")
        for name in atypes_ordered:
            _, sig, eps = LJ[name]
            lam = 1.0  # template value; FEP input sets Li+ pairs to ${LAMBDA}
            f.write(f"  {atype_idx[name]}  {eps:.5f}  {sig:.4f}  {lam:.4f}  # {name}\n")
        f.write("\n")

        f.write("Bond Coeffs  # harmonic: K r0\n\n")
        for key, (idx, (K, r0)) in btypes.items():
            f.write(f"  {idx}  {K:.1f}  {r0:.3f}  # {key[0]}-{key[1]}\n")
        f.write("\n")

        f.write("Angle Coeffs  # harmonic: K theta0\n\n")
        for key, (idx, (K, th)) in atypes_ang.items():
            f.write(f"  {idx}  {K:.1f}  {th:.1f}  # {key[0]}-{key[1]}-{key[2]}\n")
        f.write("\n")

        f.write("Dihedral Coeffs  # opls: c1 c2 c3 c4\n\n")
        for key, (idx, (c1,c2,c3,c4)) in dtypes.items():
            f.write(f"  {idx}  {c1:.3f}  {c2:.3f}  {c3:.3f}  {c4:.3f}"
                    f"  # center: {key[0]}-{key[1]}\n")
        f.write("\n")

        f.write("Atoms  # full: id mol type charge x y z\n\n")
        for gid, (mol_id, atype, q, x, y, z) in enumerate(all_atoms, start=1):
            f.write(f"  {gid}  {mol_id}  {atype}  {q:.4f}"
                    f"  {x:.6f}  {y:.6f}  {z:.6f}\n")
        f.write("\n")

        f.write("Bonds\n\n")
        for bid, (btype, i, j) in enumerate(all_bonds, start=1):
            f.write(f"  {bid}  {btype}  {i}  {j}\n")
        f.write("\n")

        f.write("Angles\n\n")
        for aid, (atype_a, i, j, k) in enumerate(all_angles, start=1):
            f.write(f"  {aid}  {atype_a}  {i}  {j}  {k}\n")
        f.write("\n")

        f.write("Dihedrals\n\n")
        for did, (dtype, i, j, k, l) in enumerate(all_dihe, start=1):
            f.write(f"  {did}  {dtype}  {i}  {j}  {k}  {l}\n")
        f.write("\n")

    n_solvent_atoms = n_mol * n_at
    print(f"  {out_path.name}: {n_atoms_total} atoms ({n_mol}×{n_at}+1), "
          f"box={L_box:.1f} Ang, "
          f"{n_bonds_total} bonds, {n_angl_total} angles, {n_dihe_total} dihedrals")


# ─── Main ────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    DATA = pathlib.Path('tools/fep/data')

    # 60 EC + 1 Li+  (target L=19.8 Ang × 1.5 = 29.8 Ang initial)
    print("Building EC system ...")
    build_system(ec_template,  60, DATA/'ec_li_aa.data',  'EC',  L_scale=1.5)

    # 50 PC + 1 Li+  (target L=19.1 Ang × 1.5 = 28.6 Ang initial)
    print("Building PC system ...")
    build_system(pc_template,  50, DATA/'pc_li_aa.data',  'PC',  L_scale=1.5)

    # 50 DME + 1 Li+ (target L=19.2 Ang × 1.5 = 28.8 Ang initial)
    print("Building DME system ...")
    build_system(dme_template, 50, DATA/'dme_li_aa.data', 'DME', L_scale=1.5)

    print("Done.")
