output_file = "../output/customplot.pdf"

# create plot
set terminal pdfcairo
set output output_file
set ylabel "Solution"
set xlabel "x"
plot "../output/solution_Gal.gpl" using 1:2 title "no viscosity" with lines,\
     "../output/solution_low.gpl" using 1:2 title "low-order viscosity" with lines,\
     "../output/solution_EV.gpl" using 1:2 title "entropy viscosity" with lines
