input_file = "../output/energy-001.csv"
output_file = "../output/energy.pdf"
y_label = "Energy"

set datafile separator ","

# create plot
set terminal postscript enhanced color
set output "| ps2pdf - ".output_file
set ylabel y_label
set xlabel "x"
plot input_file using 1:2 notitle with linesp
reset
set output
set terminal pop
