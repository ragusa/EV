set terminal postscript enhanced color
set output '| ps2pdf - compare_viscosity.pdf'
set ylabel "Viscosity"
set xlabel "x"
set title "Comparison of Entropy Viscosity"
set key bottom right
set style line 1 linetype 1 linecolor 1 pointtype 1
set style line 2 linetype 1 linecolor 3 pointtype 2
plot "entropy_viscosity_dealii.txt" using 1:2 title "deal.ii" with linesp ls 1,\
     "entropy_viscosity_MATLAB.txt" using 1:2 title "MATLAB"   with linesp ls 2
reset
set output
set terminal pop
