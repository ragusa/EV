set terminal postscript enhanced color
set output '| ps2pdf - compare_solutions.pdf'
set ylabel "Solution"
set xlabel "x"
set title "Comparison of Solutions"
set key bottom right
set style line 1 linetype 1 linecolor 1
set style line 2 linetype 1 linecolor 2
set style line 3 linetype 1 linecolor 3
set style line 4 linetype 1 linecolor 4
plot "solutions_dealii.txt" using 1:2 title "deal.ii new" with linesp ls 1,\
     "solutions_dealii.txt" using 1:3 title "deal.ii old" with linesp ls 2,\
     "solutions_MATLAB.txt" using 1:2 title "MATLAB new" with linesp ls 3,\
     "solutions_MATLAB.txt" using 1:3 title "MATLAB old" with linesp ls 4
reset
set output
set terminal pop
