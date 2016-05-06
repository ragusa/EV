# note: this script requires Gnuplot version 4.6 or higher
# usage: gnuplot -e 'problem_name=<problem name>;
#           timeintegrator=<time integrator string>'
#           viscosity.gp

filebase = "viscosity"

# list of possible input files to plot and their corresponding titles
outdir = "../output/".problem_name."/"
input_file = outdir."viscosity_EV_" .timeintegrator.".gpl"

set terminal pdfcairo
output_file = outdir."viscosity_".timeintegrator.".pdf"
set output output_file
set ylabel "Viscosity"
set xlabel "x"
set logscale y
set key top right
set style fill noborder

plot input_file using 1:2 with filledcurves y1=0 fs transparent solid 0.75\
  linetype 1 linecolor rgb "forest-green" title "Low-Order",\
     input_file using 1:3 with filledcurves y1=0 fs transparent solid 0.5\
  linetype 1 linecolor rgb "gold" title "Entropy",\
     input_file using 1:4 with filledcurves y1=0 fs transparent solid 0.25\
  linetype 1 linecolor rgb "red" title "High-Order"
