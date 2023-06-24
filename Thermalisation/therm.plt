set terminal pdfcairo size 8cm,4cm

set grid

E = 1

array y1[6]
y1[1] = 0.73
y1[2] = 0.192
y1[3] = 0.04
y1[4] = 0.005
y1[5] = -0.002
y1[6] = -0.002


array y2[6]
y2[1] = 0.74
y2[2] = 0.202
y2[3] = 0.05
y2[4] = 0.015
y2[5] = 0.008
y2[6] = 0.008


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


array c[8]
c[1] = 1000
c[2] = 2000
c[3] = 5000
c[4] = 10000
c[5] = 20000
c[6] = 30000
c[7] = 40000
c[8] = 50000

M = 1.7

set xlabel "MC Time"

do for [j=1:6:1] {
    file_cold = sprintf("24x24x24x24_b3.09_WL_%d_%d_cold.dat", loopr[j], loopt[j])
    file_hot = sprintf("24x24x24x24_b3.09_WL_%d_%d_hot.dat", loopr[j], loopt[j])

    if (j==6) {
        stats [5000:10000][*:*] file_cold using ($0):1 name "A" nooutput
        max = A_mean_y + A_stddev_y*M
        min = A_mean_y - A_stddev_y*M
    }
    else {
        stats [30000:50000][*:*] file_cold using ($0):1 name "A" nooutput
        max = A_mean_y + A_stddev_y*M
        min = A_mean_y - A_stddev_y*M
    }

    print A_mean_y, A_mean_err_y

    f(x) = A_mean_y
    
    #set title sprintf("WL(%d,%d)", loopr[j], loopt[j])
    set ytics y1[j], 0.002, y2[j]
    set yrange [y1[j]:y2[j]]
    
    do for [k=1:8:1] {
        set output sprintf("../Latex/24x24x24x24_b3.09_WL_%d_%d_%d.pdf", loopr[j], loopt[j], c[k])
        set xrange[0:c[k]]
        plot    file_hot u 0:1 w l lt rgb "#FA9000" title "hot start",\
                file_cold u 0:1 w l lt rgb "#0090FA" title "cold start",\
                f(x) lt rgb "#000000" title "",\
                min w filledcurves y1=max lc rgb "#88444444" title ""\
                
    }
}

file_cold = "24x24x24x24_b3.09_WL_10_16_cold_NAPE_100.dat"
file_hot = "24x24x24x24_b3.09_WL_10_16_hot_NAPE_100.dat"

stats [5000:10000][*:*] file_cold using ($0):1 name "A" nooutput
max = A_mean_y + A_stddev_y*M
min = A_mean_y - A_stddev_y*M

print A_mean_y, A_mean_err_y
f(x) = A_mean_y

set ytics y1[6], 0.002, y2[6]
set yrange [y1[6]:y2[6]]
    
do for [k=1:8:1] {
    set output sprintf("../Latex/24x24x24x24_b3.09_WL_10_16_%d_NAPE_100.pdf", c[k])
    set xrange[0:c[k]]
    plot    file_hot u 0:1 w l lt rgb "#FA9000" title "hot start",\
            file_cold u 0:1 w l lt rgb "#0090FA" title "cold start",\
            f(x) lt rgb "#000000" title "",\
            min w filledcurves y1=max lc rgb "#88444444" title ""\
                
}








#show variables A_
