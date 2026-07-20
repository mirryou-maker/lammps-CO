# Verification decks for styles the force-style suite does not reach

LAMMPS's `unittest/force-styles` regression suite covers most of the pair styles
modified in Step A-3, but eleven have no test in that suite and no input deck in
`examples/`. These decks close that gap: each exercises one such style, and
`compare_decks.sh` runs it under an unmodified build and an A-3 build of the same
LAMMPS commit, then compares the thermo output verbatim.

## Usage

```sh
bash compare_decks.sh decks.tsv     # <style> <TAB> <deck path> [<TAB> workdir]
```

`strip_timing.sed` removes the lines that legitimately differ between two runs of
the same physics (wall-clock numbers, the per-section timing table, host banners);
everything else must match character for character.

## Coverage

| style | note |
|---|---|
| `brownian`, `brownian/poly` | `atom_style sphere`; `brownian/poly` needs `newton off` |
| `dsmc` | |
| `eam/apip` | `atom_style apip`; reads a standard EAM setfl file |
| `multi/lucy/rx`, `table/rx` | adapted from `examples/PACKAGES/dpd-react/dpdrx-shardlow`, which supplies the `fix rx` kinetics; `rxn.table` is generated for `table/rx` |
| `polymorphic` | `CuTa_eam.poly` from `potentials/` |
| `sph/heatconduction`, `sph/lj` | `atom_style sph` |

`sph/idealgas` and `sph/taitwater/morris` are covered by the decks already in
`examples/PACKAGES/sph/` and need nothing here.
