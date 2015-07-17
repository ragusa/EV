set terminal postscript enhanced color
set output '| ps2pdf - FCT_comparison_ref8.pdf'

# begin multiplot, used here to overlay one plot on another
set multiplot

# first (base) plot
set size 1,1
set origin 0,0
set ylabel "Solution"
set xlabel "x"
set key top right
set style line 1 linetype 1 linecolor -1 pointtype -1
set style line 2 linetype 1 linecolor 1 pointtype -1
set style line 3 linetype 1 linecolor 3 pointtype -1
# draw box
set object 1 rect from 0.3,0.3 to 0.35,0.35 fc rgb "black"\
  fs pattern 0 bo 9
# draw arrows
set arrow 1 from 0.3,0.35 to 0.65,0.415 lt 2 lc 9 nohead
set arrow 2 from 0.3,0.3 to 0.65,0.21 lt 2 lc 9 nohead
plot "exact.gpl" using 1:2 title "Exact" with linesp ls 1,\
     "GalFCT_ref8.gpl" using 1:2 title "Gal-FCT" with linesp ls 2,\
     "EVFCT_ref8.gpl" using 1:2 title "EV-FCT" with linesp ls 3
set noarrow

# second plot
set size 0.3,0.45
set origin 0.65,0.4
set xrange [0.3:0.35]
set xtics 0.05
set mxtics 5
unset ytics
#set ytics off
#set ytics format ""
set format y ""
set ylabel ""
set xlabel ""
set style line 1 linetype 1 linecolor -1 pointtype -1
set style line 2 linetype 1 linecolor 1 pointtype -1
set style line 3 linetype 1 linecolor 3 pointtype -1
plot "exact.gpl" using 1:2 title "" with linesp ls 1,\
     "GalFCT_ref8.gpl" using 1:2 title "" with linesp ls 2,\
     "EVFCT_ref8.gpl" using 1:2 title "" with linesp ls 3

# end multiplot
unset multiplot

reset
set output
set terminal pop
