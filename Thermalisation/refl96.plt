set terminal pdfcairo size 8cm,4cm

set grid

E = 1

array y1[3]
y1[1] = 0.044
y1[2] = 0.009
y1[3] = 0.0015


array y2[3]
y2[1] = 0.045
y2[2] = 0.01
y2[3] = 0.0025


array loopr[3]
loopr[1] = 5
loopr[2] = 7
loopr[3] = 9

array loopt[3]
loopt[1] = 5
loopt[2] = 7
loopt[3] = 9

set xlabel "MC Time"

set xrange[0:300]
#set xtics 0,100,300

f(x) = c

M = 0.12

do for [j=1:3:1] {
    set output sprintf("../Latex/96x96x96x96_b3.09_WL_%d_%d.pdf", loopr[j], loopt[j])
    
    file_1 = sprintf("24x24x24x24_b3.09_WL_%d_%d_cold.dat", loopr[j], loopt[j])
    file_2 = sprintf("96x96x96x96_b3.09_WL_%d_%d.dat", loopr[j], loopt[j])

    stats [30000:50000][*:*] file_1 using ($0):1 name "B"

    c = B_mean_y
    max = B_mean_y + B_stddev_y*M
    min = B_mean_y - B_stddev_y*M

    set ytics y1[j], 0.001, y2[j]
    set yrange [y1[j]:y2[j]]
    
    plot    min w filledcurves y1=max lc rgb "#99444444" title "",\
            file_2 u ($0):1 w l lt rgb "#9900CC" title "",\
            f(x) lt rgb "#000000" title ""
}



