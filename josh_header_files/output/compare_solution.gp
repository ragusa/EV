output_file = "compare_solution.pdf"

# create plot
set terminal postscript enhanced color
set output "| ps2pdf - ".output_file
set ylabel "Solution"
set xlabel "x"
plot "solution_none.gpl" using 1:2 title "none" with lines,\
     "solution_const.gpl" using 1:2 title "const" with lines,\
     "solution_first.gpl" using 1:2 title "first order" with lines,\
     "solution_ent.gpl" using 1:2 title "entropy" with lines,\
     "solution_jumps.gpl" using 1:2 title "entropy with jumps" with lines
reset
set output
set terminal pop
