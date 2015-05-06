set terminal postscript enhanced color
set output '| ps2pdf - speedup.pdf'
set ylabel "Speedup"
set xlabel "Number of Processes"
set key top left
set style line 1 linetype 1 linecolor 0
set style line 2 linetype 1 linecolor 1 pointtype 1
set style line 3 linetype 1 linecolor 2 pointtype 1
set style line 4 linetype 1 linecolor 3 pointtype 1
set style line 5 linetype 1 linecolor 4 pointtype 1
set style line 6 linetype 1 linecolor 5 pointtype 1
set xrange [0:25]
linear(x) = x
plot linear(x) title "Linear Reference" with lines ls 1,\
   "speedup_initialize.gpl" title "Initialize" with errorbars ls 2,\
   "speedup_setup.gpl"      title "Setup"      with errorbars ls 3,\
   "speedup_assemble.gpl"   title "Assemble"   with errorbars ls 4,\
   "speedup_solve.gpl"      title "Solve"      with errorbars ls 5,\
   "speedup_output.gpl"     title "Output"     with errorbars ls 6,\
   "speedup_initialize.gpl" using 1:2 notitle with lines ls 2,\
   "speedup_setup.gpl"      using 1:2 notitle with lines ls 3,\
   "speedup_assemble.gpl"   using 1:2 notitle with lines ls 4,\
   "speedup_solve.gpl"      using 1:2 notitle with lines ls 5,\
   "speedup_output.gpl"     using 1:2 notitle with lines ls 6
