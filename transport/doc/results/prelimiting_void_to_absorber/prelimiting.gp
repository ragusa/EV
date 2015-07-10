set terminal postscript enhanced color
set output '| ps2pdf - prelimiting.pdf'
set ylabel "Solution"
set xlabel "x"
set key bottom left
set style line 1 linetype 1 linecolor 0
set style line 2 linetype 1 linecolor 3 pointtype 1
set style line 3 linetype 1 linecolor 1 pointtype 6
set style line 4 linetype 2 linecolor 2 pointtype 2
set style line 5 linetype 2 linecolor 4 pointtype 4
plot "uexact.txt" using 1:2 title "Exact" with lines ls 1,\
     "uFCT_Gal_FE_prelimit.txt" using 1:2 title "FE, prelimiting" with linesp ls 3,\
     "uFCT_Gal_FE_noprelimit.txt" using 1:2 title "FE, no prelimiting" with linesp ls 2,\
     "uFCT_Gal_BE.txt" using 1:2 title "BE, no prelimiting" with linesp ls 4,\
     "uFCT_Gal_SSPRK33.txt" using 1:2 title "SSPRK33, no prelimiting" with linesp ls 5
reset
set output
set terminal pop
