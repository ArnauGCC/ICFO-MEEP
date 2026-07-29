#!/bin/bash

SCRIPT="simulation.py"

for NP in 8 14 16
do
    echo "===== np=$NP ====="

    for i in {1..3}
    do
        /usr/bin/time -f "Run $i: %e s" \
            mpirun -np $NP python $SCRIPT >/dev/null
    done
done