input_file = "momentum-001.csv"
output_file = "momentum.pdf"
y_label = "Momentum"

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
