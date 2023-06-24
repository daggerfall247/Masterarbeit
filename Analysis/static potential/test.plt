set terminal pdfcairo size 8cm,5cm

a = 0.012899
u = 5.067731163
o = 2.7225

a_carolin_c = 0.026
offset_c = 2.7225

set output "test.pdf"

set ylabel "V(R=12) [GeV]"

set xlabel "1/Vol [1/fm^4]"
#set xlabel "L [fm]"

set xrange [*:*]

set key bottom right
set pointsize 0.3

V = 6.1
m = 1.61

f(x) = V + A/x*exp(-sqrt(3)*m/u*x/2)

g(x) = V + x/(2*y)*M2

fit g(x) "R12.dat" using ($2/(a**4)):($3/(a*u)):($4/(a*u)) via V, y, M2
plot    g(x), "R12.dat" using ($2/(a**4)):($3/(a*u)):($4/(a*u)) w yerrorbars title ""

#fit f(x) "R12.dat" using (1/($2**(0.25))*a):($3/(a*u)) via A,m,V
#plot    f(x), "R12.dat" using (1/($2**(0.25))*a):($3/(a*u)):($4/(a*u)) w yerrorbars title ""
