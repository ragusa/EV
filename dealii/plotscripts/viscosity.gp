output_file = "../output/viscosity.pdf"
set datafile separator ","

# create plot
set terminal postscript enhanced color
set output "| ps2pdf - ".output_file
set ylabel "Viscosity"
set xlabel "x"
set logscale y
plot "../output/viscosity.csv" using 1:2 title "used" with lines,\
"../output/first_order_viscosity.csv" using 1:2 title "first order" with lines,\
"../output/entropy_viscosity.csv" using 1:2 title "entropy" with lines,\
"../output/entropy_viscosity_with_jumps.csv" using 1:2 title "entropy with jumps" with lines
reset
set output
set terminal pop
