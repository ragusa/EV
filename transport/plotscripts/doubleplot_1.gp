set terminal postscript enhanced color
set output '| ps2pdf - ../plots/test_problem_1.pdf'
set multiplot layout 2, 1 title "Test Problem 1, t = 0.5"
set ylabel "Angular Flux"
set key bottom right
set bmargin 0
set format x ""
set format y "%6g"
plot "../output/exact_solution_1.gpl" using 1:2 title "Exact" with lines,\
     "../output/solution_low_order_1.gpl" using 1:2 title "Low-Order" with linesp,\
     "../output/solution_high_order_1.gpl" using 1:2 title "High-Order" with linesp
set title ""
set ylabel "Viscosity"
set xlabel "x"
#set logscale y
set bmargin
set tmargin 0
set format x "%g"
set format y "%6g"
#set yrange [0:60]
#set ytics 0,10,50
set key top right
plot "../output/viscosity_1.gpl" using 1:2 title "Low-Order" with linesp,\
     "../output/viscosity_1.gpl" using 1:3 title "Entropy"   with linesp
unset multiplot
