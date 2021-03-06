# note: this script requires Gnuplot version 4.6 or higher
# usage: gnuplot -e 'problem_name=<problem name>;
#           timeintegrator=<time integrator string>'
#           viscosity.gp

filebase = "viscosity"

# list of possible input files to plot and their corresponding titles
outdir = "../output/".problem_name."/"
input_file = outdir."viscosity_EV_" .timeintegrator.".gpl"

set terminal postscript enhanced color
output_file = outdir."viscosity_".timeintegrator.".pdf"
set output '| ps2pdf - '.output_file
set ylabel "Viscosity"
set xlabel "x"
set logscale y
set key top right

plot input_file using 1:2 with lines\
  linetype 2 linecolor rgb "green" title "Low-Order Viscosity",\
     input_file using 1:3 with lines\
  linetype 1 linecolor rgb "blue"  title "Entropy Viscosity"
