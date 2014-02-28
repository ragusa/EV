input_file = "../output/burgers-001.gpl"
exact_file = "../output/burgers_exact-001.gpl"
output_file = "../output/burgers.pdf"
y_label = "Velocity"

# create plot
set terminal postscript enhanced color
set output "| ps2pdf - ".output_file
set ylabel y_label
set xlabel "x"
#plot input_file using 1:2 title "FEM Solution" with linesp
plot input_file using 1:2 title "FEM Solution" with linesp,\
     exact_file using 1:2 title "Exact Solution" with linesp
reset
set output
set terminal pop
