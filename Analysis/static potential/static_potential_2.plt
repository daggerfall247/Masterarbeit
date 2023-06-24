set terminal pdfcairo size 14cm,8cm

set output "../../Latex/V_lat.pdf"

set key bottom right

set pointsize 0.3

set multiplot

set ylabel "aV(R/a)"
set xlabel "R/a"

set xrange [0.5:12.5]
set yrange [0.22:0.41]

plot    "V_24.dat" u ($1):($3):($4) w yerrorbars title "A",\
        "V_48.dat" u ($1):($3):($4) w yerrorbars title "B",\
        "V_96.dat" u ($1):($3):($4) w yerrorbars title "C",\

unset key
set origin 0.8,0.45
set size 0.12,0.4
clear

set yrange [0.389:0.397]
set xrange [10.8:11.2]
unset ylabel
unset xlabel
unset ytics
set xtics 11

plot    "V_24.dat" u ($1):($3):($4) w yerrorbars title "A",\
        "V_48.dat" u ($1):($3):($4) w yerrorbars title "B",\
        "V_96.dat" u ($1):($3):($4) w yerrorbars title "C",\

unset multiplot
