set terminal pdfcairo size 14cm,8cm

set output "../../Latex/V.pdf"

set multiplot

a = 0.012899
u = 5.067731163

a_carolin_a = 0.078
a_carolin_b = 0.041
a_carolin_c = 0.026

offset_a = 2.15
offset_b = 1.11
offset_c = 0

offset = 2.726

set key bottom right



set ylabel "V(R) [GeV]"
set xlabel "R [fm]"

#set xrange [*:0.17]


set pointsize 0.3


f(x) = A/x + B*x + C

fit f(x) "carolin_c.dat" u ($1*a_carolin_c):($2/(a_carolin_c*u) + offset_c + offset):($3/(a_carolin_c*u)) via A,B,C


plot    "V_24.dat" u ($1*a):($3/(a*u)):($4/(a*u)) w yerrorbars title "A",\
        "V_48.dat" u ($1*a):($3/(a*u)):($4/(a*u)) w yerrorbars title "B",\
        "V_96.dat" u ($1*a):($3/(a*u)):($4/(a*u)) w yerrorbars title "C",\
        "carolin_a.dat" u ($1*a_carolin_a):($2/(a_carolin_a*u) + offset_a + offset):($3/(a_carolin_a*u)) w yerrorbars title "ref. a",\
        "carolin_b.dat" u ($1*a_carolin_b):($2/(a_carolin_b*u) + offset_b + offset):($3/(a_carolin_b*u)) w yerrorbars title "ref. b",\
        "carolin_c.dat" u ($1*a_carolin_c):($2/(a_carolin_c*u) + offset_c + offset):($3/(a_carolin_c*u)) w yerrorbars title "ref. c",\
        f(x) w l lc rgb "#BB222222" title "fit V(R) on ref c"

set origin 0.29,0.21
set size 0.5,0.55
clear

unset key
unset ylabel
unset xlabel
set xrange [0.1:0.16]

plot    "V_24.dat" u ($1*a):($3/(a*u)):($4/(a*u)) w yerrorbars title "A",\
        "V_48.dat" u ($1*a):($3/(a*u)):($4/(a*u)) w yerrorbars title "B",\
        "V_96.dat" u ($1*a):($3/(a*u)):($4/(a*u)) w yerrorbars title "C",\
        "carolin_a.dat" u ($1*a_carolin_a):($2/(a_carolin_a*u) + offset_a + offset):($3/(a_carolin_a*u)) w yerrorbars title "ref. a",\
        "carolin_b.dat" u ($1*a_carolin_b):($2/(a_carolin_b*u) + offset_b + offset):($3/(a_carolin_b*u)) w yerrorbars title "ref. b",\
        "carolin_c.dat" u ($1*a_carolin_c):($2/(a_carolin_c*u) + offset_c + offset):($3/(a_carolin_c*u)) w yerrorbars title "ref. c",\
        f(x) w l lc rgb "#BB222222" title "fit V(R) on ref c"

unset multiplot


