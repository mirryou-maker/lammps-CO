# Match hybrid pair styles by prefix so an active suffix does not break sub-style checks

## Summary

`hybrid`, `hybrid/overlay`, and `hybrid/scaled` are each registered under two
names — the plain one and a `/omp` alias bound to the same class:

```
src/pair_hybrid.h:16          PairStyle(hybrid,PairHybrid);
src/pair_hybrid.h:17          PairStyle(hybrid/omp,PairHybrid);
src/pair_hybrid_overlay.h:16  PairStyle(hybrid/overlay,PairHybridOverlay);
src/pair_hybrid_overlay.h:17  PairStyle(hybrid/overlay/omp,PairHybridOverlay);
src/pair_hybrid_scaled.h:16   PairStyle(hybrid/scaled,PairHybridScaled);
src/pair_hybrid_scaled.h:17   PairStyle(hybrid/scaled/omp,PairHybridScaled);
```

With a suffix active, `Force::new_pair()` finds the aliased name, so `sflag` is
1 and `Force::store_style()` records `force->pair_style` as
`"hybrid/overlay/omp"`. Three call sites compare that field for exact equality
and therefore misbehave under `-sf omp`.

## What breaks

**1. `pair_style hybrid/overlay kolmogorov/crespi/z` cannot run under `-sf omp`.**

`PairKolmogorovCrespiZ::settings()` requires the enclosing style to be
`hybrid/overlay` and aborts when the comparison fails — which it does under a
suffix, even though `hybrid/overlay` is exactly what is in use:

```
ERROR: ERROR: requires hybrid/overlay pair_style (src/INTERLAYER/pair_kolmogorov_crespi_z.cpp:207)
Last input line: pair_style hybrid/overlay kolmogorov/crespi/z 16.0
```

**2. `compute fep` and `fix adapt/fep` silently skip a validation.**

Both use the comparison to decide whether to check that the requested type range
is served by the named sub-style. Under a suffix the comparison is false, so the
check is skipped and an invalid range is no longer reported. This is a lost
diagnostic rather than a crash, which is why it has gone unnoticed.

## Why the test suite does not catch it

`unittest/force-styles` skips its OMP branch when no `/omp` variant of the
sub-style exists, and none exists for `kolmogorov/crespi/z`, so the path is
never exercised. I ran into this while adding OMP variants for styles that
currently lack them — adding the variant is what made the branch reachable.

## The fix

Replace the exact comparison with a prefix match, which is what LAMMPS already
uses elsewhere for the same question. `Force::pair_match()` gates its own
`dynamic_cast<PairHybrid *>` on exactly this test:

```cpp
else if (utils::strmatch(pair_style, "^hybrid")) {
  auto *hybrid = dynamic_cast<PairHybrid *>(pair);
```

so the patch brings these three sites in line with that idiom rather than
introducing a new one.

**One deliberate widening.** In `compute_fep` and `fix_adapt_fep` the original
test named `hybrid` and `hybrid/overlay` but not `hybrid/scaled` or
`hybrid/molecular`, so the validation was already being skipped for those. Using
`^hybrid` covers them. Every style matched is a `PairHybrid` subclass, so the
`dynamic_cast` remains valid. If you would rather keep the fix strictly
behaviour-preserving I am happy to restrict the match to the original two names
plus their suffixed forms — please say which you prefer.

## Reproducer

```
units metal
atom_style atomic
boundary p p p
region box block 0 10 0 10 0 10
create_box 2 box
pair_style hybrid/overlay kolmogorov/crespi/z 16.0
```

```
$ lmp -in in.kc            # runs
$ lmp -in in.kc -sf omp    # ERROR: requires hybrid/overlay pair_style
```

After the patch both invocations run.

## Testing

Built from `91d4111a9` with `PKG_OPENMP KSPACE EXTRA-PAIR MANYBODY MOLECULE FEP
INTERLAYER` and `ENABLE_TESTING=on`, GCC 11.2, OpenMPI 4.1.4:

- compiles with no new warnings or errors
- `ctest -R PairStyle` — **325/325 pass**, unchanged from the same build without
  the patch
- the reproducer above runs under `-sf omp`

## Notes

No documentation change seems warranted: the fix restores the behaviour the
existing documentation already describes. Happy to add a note if you disagree.
