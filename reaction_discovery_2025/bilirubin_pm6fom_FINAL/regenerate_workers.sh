#!/usr/bin/bash
#
# regenerate_workers.sh - create secondary1..bbb from the oldSecondary template.
#
# Each worker carries its own index in the claim routine, so that worker i ignores lattice points
# below index i instead of every worker contending for the same early points:
#
#     secondary<i>/s1.py    if current == 0 or current < <i>:     for i = 1 .. bbb-1
#                           if current == 0:                      for i = bbb
#
# loadFoundlings.sh, in the same directory, performs the clone without that edit.
#
# Usage:  ./regenerate_workers.sh [bbb]      (default 400; must match bbb in p.f90)

set -eu
bbb=${1:-400}
here=$(cd "$(dirname "$0")" && pwd)
cd "$here"

[ -d oldSecondary ] || { echo "error: no oldSecondary/ template in $here" >&2; exit 1; }
grep -q "if current == 0:" oldSecondary/s1.py || {
    echo "error: oldSecondary/s1.py does not contain the unstaggered template line" >&2; exit 1; }

rm -rf secondary[0-9]*
for i in $(seq 1 "$bbb"); do
    cp -R oldSecondary "secondary$i"
    if [ "$i" -lt "$bbb" ]; then
        sed -i "s/^\( *\)if current == 0:/\1if current == 0 or current < $i:/" "secondary$i/s1.py"
    fi
done

echo "regenerated secondary1..secondary$bbb from oldSecondary"
echo "  stagger applied to secondary1..secondary$((bbb - 1)); secondary$bbb left unstaggered"
echo "next: sbatch each worker (see batch30.sh)"
