set terminal pdfcairo size 8cm,4cm

set grid

E = 1

array y1[6]
y1[1] = 0.73
y1[2] = 0.194
y1[3] = 0.043
y1[4] = 0.008
y1[5] = 0.0
y1[6] = -0.002


array y2[6]
y2[1] = 0.734
y2[2] = 0.198
y2[3] = 0.047
y2[4] = 0.012
y2[5] = 0.004
y2[6] = 0.002


array loopr[6]
loopr[1] = 1
loopr[2] = 3
loopr[3] = 5
loopr[4] = 7
loopr[5] = 9
loopr[6] = 10

array loopt[6]
loopt[1] = 1
loopt[2] = 3
loopt[3] = 5
loopt[4] = 7
loopt[5] = 9
loopt[6] = 16

unset ytics
set xlabel "MC Time"

set xrange[0:300]
#set xtics 0,100,300

f(x) = c

M = 1.7

do for [j=1:6:1] {
    set output sprintf("../Latex/48x48x48x48_b3.09_WL_%d_%d.pdf", loopr[j], loopt[j])
    
    file_1 = sprintf("24x24x24x24_b3.09_WL_%d_%d_cold.dat", loopr[j], loopt[j])
    file_2 = sprintf("48x48x48x48_b3.09_WL_%d_%d.dat", loopr[j], loopt[j])

    stats [300:1000][*:*] file_2 using ($0):1 name "A"

    if(j==6){
        stats [5000:10000][*:*] file_1 using ($0):1 name "B"
    }
    else {
        stats [30000:50000][*:*] file_1 using ($0):1 name "B"
    }

    max = B_mean_y + A_stddev_y*M
    min = B_mean_y - A_stddev_y*M
    c = B_mean_y

    set ytics y1[j], 0.001, y2[j]
    set yrange [y1[j]:y2[j]]
    
    plot    min w filledcurves y1=max lc rgb "#99444444" title "",\
            file_2 u ($0):1 w l lt rgb "#9900CC" title "",\
            f(x) lt rgb "#000000" title ""
}



