set terminal pdfcairo size 7cm,4cm

unset key
#file = "24x24x24x24_SU2_b3.09_qtop_exception_1.meas"
file = "48x48x48x48_SU2_b3.09_qtop_exception_2.meas"
set yrange [-5:5]

set mytics 2
set mxtics 5

set pointsize 0.3

set xlabel "N_{APE}"
set ylabel "Q_{TOP}"

#set output "../../Latex/qtop_48_exception_1.pdf"
set output "../../Latex/qtop_48_exception_2.pdf"

plot    "<(sed -n '/^QTOP 45/p' 48x48x48x48_SU2_b3.09_qtop_exception_2.meas)" u 3:5 w l,\
        "<(sed -n '/^QTOP 46/p' 48x48x48x48_SU2_b3.09_qtop_exception_2.meas)" u 3:5 w l,\
        "<(sed -n '/^QTOP 47/p' 48x48x48x48_SU2_b3.09_qtop_exception_2.meas)" u 3:5 w l,\
        "<(sed -n '/^QTOP 48/p' 48x48x48x48_SU2_b3.09_qtop_exception_2.meas)" u 3:5 w l,\
        "<(sed -n '/^QTOP 49/p' 48x48x48x48_SU2_b3.09_qtop_exception_2.meas)" u 3:5 w l,\
