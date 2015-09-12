output_file = "../output/customplot.pdf"

# create plot
set terminal postscript enhanced color
set output "| ps2pdf - ".output_file
set ylabel "Solution"
set xlabel "x"
plot "../output/solution_none.gpl" using 1:2 title "none" with lines,\
     "../output/solution_const.gpl" using 1:2 title "const" with lines,\
     "../output/solution_first.gpl" using 1:2 title "first order" with lines,\
     "../output/solution_ent.gpl" using 1:2 title "entropy" with lines,\
     "../output/solution_jumps.gpl" using 1:2 title "entropy with jumps" with lines
