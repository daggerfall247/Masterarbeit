#!/bin/bash

for r in {1..12}
do
    truncate -s 0 "R$r.dat"
    for i in 24 48 96
    do
        sed -n "${r}p" V_$i.dat >> "R$r.dat"
    done
done
