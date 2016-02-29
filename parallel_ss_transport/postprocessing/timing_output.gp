set terminal postscript enhanced color
set output '| ps2pdf - timing_output.pdf'
set ylabel "Wallclock Time"
set xlabel "Problem Size"
set key top left
set logscale x
set logscale y
set style line 1 linetype 1 linecolor 0 pointtype 1
set style line 2 linetype 1 linecolor 1 pointtype 1
set style line 3 linetype 1 linecolor 2 pointtype 1
set style line 4 linetype 1 linecolor 3 pointtype 1
plot \
   "timing_p00_output.gpl" title "Serial" with errorbars ls 1,\
   "timing_p01_output.gpl" title "p = 1" with errorbars ls 2,\
   "timing_p08_output.gpl" title "p = 8" with errorbars ls 3,\
   "timing_p20_output.gpl" title "p = 20" with errorbars ls 4,\
   "timing_p00_output.gpl" notitle with lines ls 1,\
   "timing_p01_output.gpl" notitle with lines ls 2,\
   "timing_p08_output.gpl" notitle with lines ls 3,\
   "timing_p20_output.gpl" notitle with lines ls 4
