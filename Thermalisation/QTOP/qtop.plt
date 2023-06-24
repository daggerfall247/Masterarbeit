set terminal pdfcairo

alpha = 0.4

set output sprintf("qtop_alpha_%f.pdf", alpha)
set grid

set xlabel "N_{APE}"
set ylabel "Q_{TOP}"

set offsets graph 0, 0, 0.5, 0.5

do for [k=10000:50000:2500] {
    file = sprintf("24x24x24x24_SU2_b3.09_hot_qtop_%d.dat", k)
    
    set title sprintf("Config: %d, {/Symbol a}_{APE} = %f", k, alpha)
    plot file u 1:($3) w l title ""
}

#plot "QTOP.dat" u 1:4 w lp title ""
