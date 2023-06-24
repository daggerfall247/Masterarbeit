set terminal pdfcairo size 14cm,3cm

unset key
set grid mxtics mytics xtics ytics
 
set pointsize 0.3

set yrange [-2.4:2.4]
set ytics -5,1,5



set ylabel "Q_{TOP}"
set xlabel "# config"

test = "this is a string"

do for [t in "24 48 96"] {
    file = sprintf("qtop_%s.dat", t)
    set output sprintf("../../Latex/qtop_%s_all.pdf", t)
    stats file using ($0):4 name "A" nooutput
    if (t=="96") {
        set xrange [-0.01:400.01];
        set xtics 0,100,400;
        set mxtics 10
    } else {
        set xrange [-0.01:1000.01];
        set xtics 0,250,1000;
        set mxtics 5
    }
    plot file u ($0):4 with p title ""
}


set terminal pdfcairo size 5cm,5cm

unset grid

set xrange [-3.4:3.4]
set yrange [0:1]

set ytics 0,0.2,1
set xtics -3,1,3

unset mxtics
unset mytics

set xlabel "Q_{TOP}"
set ylabel "N(Q)/N_{conf}"

#set key top right

set style fill solid 0.7

max = 3.5
min = -3.5

binwidth=1.0
set boxwidth 0.95*binwidth
bin(x,width)=width*(floor((x-min)/width)+0.5) + min

f(x, sigma) = 1/sqrt(2*pi*sigma**2)*exp(-x**2/(2*sigma**2))

do for [t in "24 48 96"] {
    file = sprintf("qtop_%s.dat", t)
    set output sprintf("../../Latex/qtop_%s_histo.pdf", t)
    
    if (t=="96") {
        plot file using (bin($4,binwidth)):(0.0025) smooth freq with boxes title "comp. distr.", f(x, 2.1678) title "exp. distr."
    }
    if (t=="48") {
        plot file using (bin($4,binwidth)):(0.001) smooth freq with boxes title "comp. distr." , f(x, 0.5419) title "exp. distr."
    }
    if (t=="24") {
        plot file using (bin($4,binwidth)):(0.001) smooth freq with boxes title "comp. distr." , f(x, 0.1356) title "exp. distr."
    }
}
