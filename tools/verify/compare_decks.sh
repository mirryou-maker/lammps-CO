#!/bin/bash
# Run each deck under the unmodified and the A-3 build and compare the thermo
# output verbatim. Styles reached by neither the force-style suite nor an
# examples/ deck are covered here.
#
# Usage: compare_decks.sh <deck-list-file>
#   each line: <style> <tab> <absolute path to input deck> [<tab> workdir]
set -u
. /etc/profile.d/modules.sh
module load CUDA/12.4 DEVTOOLSET/11 OPENMPI/4.1.4.GCC8.5 >/dev/null 2>&1

BASE=/home/cyyou68/repos/lammps-co-base/build-most/lmp
PATCHED=/home/cyyou68/repos/lammps-co-patched/build-most/lmp
OUT=/scratch/cyyou68/co_deck_compare
rm -rf "$OUT"; mkdir -p "$OUT"

printf "%-26s %-10s %s\n" "style" "result" "note"
printf -- "-%.0s" {1..70}; echo

pass=0; fail=0; err=0
while IFS=$'\t' read -r style deck workdir; do
    [ -z "${style:-}" ] && continue
    [ "${style:0:1}" = "#" ] && continue
    wd="${workdir:-$(dirname "$deck")}"
    tag=$(echo "$style" | tr '/' '_')

    for b in base patched; do
        exe=$BASE; [ "$b" = patched ] && exe=$PATCHED
        ( cd "$wd" && "$exe" -in "$deck" -log "$OUT/${tag}_$b.log" -screen none ) \
            > "$OUT/${tag}_$b.out" 2>&1
        echo "$?" > "$OUT/${tag}_$b.rc"
    done

    rc_b=$(cat "$OUT/${tag}_base.rc"); rc_p=$(cat "$OUT/${tag}_patched.rc")
    if [ "$rc_b" != 0 ] || [ "$rc_p" != 0 ]; then
        note=$(grep -m1 -iE "ERROR" "$OUT/${tag}_base.log" "$OUT/${tag}_patched.log" 2>/dev/null \
               | head -1 | cut -c1-40)
        printf "%-26s %-10s %s\n" "$style" "ERROR" "rc=$rc_b/$rc_p ${note:-}"
        err=$((err + 1)); continue
    fi

    # compare thermo bodies only: drop timing/host lines that legitimately differ
    for b in base patched; do
        sed -f /scratch/cyyou68/strip_timing.sed "$OUT/${tag}_$b.log" > "$OUT/${tag}_$b.thermo"
    done

    if diff -q "$OUT/${tag}_base.thermo" "$OUT/${tag}_patched.thermo" > /dev/null; then
        printf "%-26s %-10s %s\n" "$style" "IDENTICAL" ""
        pass=$((pass + 1))
    else
        n=$(diff "$OUT/${tag}_base.thermo" "$OUT/${tag}_patched.thermo" | grep -c '^[<>]')
        printf "%-26s %-10s %s\n" "$style" "DIFFER" "$n lines"
        fail=$((fail + 1))
    fi
done < "$1"

echo
echo "identical=$pass  differ=$fail  error=$err"
