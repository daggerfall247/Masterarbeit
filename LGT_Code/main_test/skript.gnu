set terminal pdfcairo color enhanced
set output "histo.pdf"

set xrange[-1.2:1.2]

set samples 5000

set key top left

alpha = 7.0
N = (pi/(2*alpha**3))**(1/3)*exp(alpha)
f(x, y, z) = 1/N * sqrt(1.0-x**2)*exp(y*x)
M = N/35.05
M2 = 4

binwidth=2/M2
V = binwidth * M

plot  f(x, alpha, N) lw 1, "improved.dat" u 1:($2/M) w l lw 0.3, "oldschool.dat" u 1:($2/M) w l lw 0.3
