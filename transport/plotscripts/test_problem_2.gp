set terminal postscript enhanced color
set output '| ps2pdf - test_problem_2.pdf'
set ylabel "Angular Flux"
set xlabel "x"
set key top right
set yrange [-0.1:1.1]
plot "analytic_solution_test_problem_2.dat" using 1:2 title "Analytic Solution" with lines,\
"solution_none_2.gpl" using 1:2 title "No Viscosity Solution" with lines,\
"solution_first_order_2.gpl" using 1:2 title "First-Order Viscosity Solution" with lines,\
"solution_entropy_2.gpl" using 1:2 title "Entropy Viscosity Solution" with lines,\
"solution_max_principle_2.gpl" using 1:2 title "Maximum Principle Viscosity Solution" with lines
reset
set output
set terminal pop
