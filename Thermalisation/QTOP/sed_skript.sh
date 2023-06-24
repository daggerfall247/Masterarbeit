#!/bin/bash

for k in {10000..50000..2500}
do
    cat "24x24x24x24_SU2_b3.090000_id1_hot_qtop.meas" | sed -n -e "s/^.*QTOP $k //p" > "24x24x24x24_SU2_b3.09_hot_qtop_$k.dat"
done

#cat "../16x16x16x16_SU2_b2.300000_id1.meas" | sed -n -e "s/^.*QTOP //p" > "QTOP.dat"
