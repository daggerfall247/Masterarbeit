set terminal pdfcairo

set grid

E = 1

set pointsize 0.3
set bars small

loop = "WL(7,7)"
set output "24x24x24x24_b3.09_WL_7_7_5000.pdf"
file_fit = "fit_5000_24x24x24x24_b3.09_WL_7_7.dat"
intervall_length = 5000





do for [x in "50000 10000 2000"] {
    set xrange[0:x]

    
    #set yrange[0.73:0.74]
    #set yrange[0.04:0.05]
    #set yrange[0:1]
    #set yrange[0.0018:0.0022]
    set title sprintf("average over an intervall of %d data points startint a t", intervall_length)
    set ylabel sprintf("<%s>_{cold} +- {/Symbol s}_{cold}", loop)
    set xlabel "t [MC time] (starting point of intervall)"
    plot file_fit u 1:2 w lines title ""
    

    set yrange[-6:6]
    set title sprintf("difference of both averages over an intervall of %d data points starting at t \ndivided by the standard error of the mean", intervall_length)
    set ylabel sprintf("(<%s>_{cold} - <%s>_{hot}) / {/Symbol s}_{cold}", loop, loop)
    set xlabel "t [MC time] (starting point of intervall)"

    plot file_fit u 1:($2-$3)/sqrt($4*$4+$5*$5) w lines title ""
    #plot file_fit u 1:($2-$3)/($5) w lines title ""
}










#file_cold = "24x24x24x24_b3.09_WL_9_9_cold.dat"
#file_hot = "24x24x24x24_b3.09_WL_9_9_hot.dat"


#do for [t=0:50000-intervall_length+1:1000]{
#   set table
#   plot file_fit u 0:($0==t?(meanCold=$2):$2), '' u 0:($0==t?(meanHot=$3):$3), '' u 0:($0==t?(errorCold=$4):$4), '' u 0:($0==t?(errorHot=$5):$5)
#   unset table
#   set title sprintf("WL(1,1), fit range = [%5d, %5d]", t, t+intervall_length)
#   f(x) = meanCold
#   g(x) = meanHot
#   plot file_cold u ($0*E):1 every E w l title "cold start", file_hot u ($0*E):1 every E w l title "hot start", f(x) lw 2 title "avg from cold start", g(x) lw 2 title "avg from hot start"
#}

#plot file_fit u 1:($2-$3):($4*1) w yerrorbars title "avg diff (cold - hot) /w cold fit error"
#plot file_fit u 1:($2-$3):($4*2) w yerrorbars title "avg diff (cold - hot) /w cold fit error"
#plot file_fit u 1:($2-$3):($4*3) w yerrorbars title "avg diff (cold - hot) /w cold fit error"

#plot file_fit u 1:($2-$3):($5*1) w yerrorbars title "avg diff (cold - hot) /w hot fit error"
#plot file_fit u 1:($2-$3):($5*2) w yerrorbars title "avg diff (cold - hot) /w hot fit error"
#plot file_fit u 1:($2-$3):($5*3) w yerrorbars title "avg diff (cold - hot) /w hot fit error"
