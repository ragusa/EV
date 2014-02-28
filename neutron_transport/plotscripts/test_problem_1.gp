set terminal postscript enhanced color
set output '| ps2pdf - test_problem_1.pdf'
set ylabel "Angular Flux"
set xlabel "x"
set key bottom right
set style line 1 linetype 1 linecolor 2
set style line 1 linetype 0 linecolor 0
plot "analytic_solution_test_problem_1.dat" using 1:2 title "Analytic Solution" with lines,\
"solution_none_1.gpl" using 1:2 title "No Viscosity Solution" with lines,\
"solution_first_order_1.gpl" using 1:2 title "First-Order Viscosity Solution" with lines,\
"solution_entropy_1.gpl" using 1:2 title "Entropy Viscosity Solution" with lines,\
"solution_max_principle_1.gpl" using 1:2 title "Maximum Principle Viscosity Solution" with lines
reset
set output
set terminal pop
