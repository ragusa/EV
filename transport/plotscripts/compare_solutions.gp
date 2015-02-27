set terminal postscript enhanced color
set output '| ps2pdf - ../plots/compare_solutions.pdf'
set ylabel "Solution"
set xlabel "x"
set title "Comparison of Solutions"
set key bottom right
set style line 1 linetype 1 linecolor 1 pointtype 1
set style line 2 linetype 1 linecolor 2 pointtype 2
plot "../output/solution_high_order_2.gpl" using 1:2 title "deal.ii" with linesp ls 1,\
     "../output/uH_EV_FE.txt" using 1:2 title "MATLAB" with linesp ls 2
reset
set output
set terminal pop
