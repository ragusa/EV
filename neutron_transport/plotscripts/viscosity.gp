set terminal postscript enhanced color
set output '| ps2pdf - viscosity.pdf'
set ylabel "Artificial Viscosity"
set xlabel "x"
set key bottom right
set logscale y
set style line 1 linetype 1 linecolor 2
set style line 2 linetype 1 linecolor 1
#set yrange [1e-5:1]
plot "viscosity.gpl" using 1:2 title "Max viscosity" with lines ls 1,\
"viscosity.gpl" using 1:3 title "Entropy viscosity" with linesp ls 2
reset
set output
set terminal pop
