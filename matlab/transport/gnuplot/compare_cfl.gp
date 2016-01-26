output_file = "../output/compare_cfl.pdf"

# create plot
set terminal postscript enhanced color
set output '| ps2pdf - '.output_file
set ylabel "Solution"
set xlabel "x"
plot "../output/uexact.txt" using 1:2 linecolor 0 pointtype 0 title "Exact" with linesp,\
     "../output/uFCT_Gal_BE_CFL0.1.txt" using 1:2 linecolor 1 pointtype 1\
       linetype 2 title "CFL=0.1" with linesp,\
     "../output/uFCT_Gal_BE_CFL0.25.txt" using 1:2 linecolor 2 pointtype 2\
       linetype 3 title "CFL=0.25" with linesp,\
     "../output/uFCT_Gal_BE_CFL0.5.txt" using 1:2 linecolor 3 pointtype 3\
       linetype 4 title "CFL=0.5" with linesp,\
     "../output/uFCT_Gal_BE_CFL1.txt" using 1:2 linecolor 4 pointtype 4\
       linetype 5 title "CFL=1.0" with linesp,\
     "../output/uFCT_Gal_BE_CFL5.txt" using 1:2 linecolor 1 pointtype 6\
       linetype 6 title "CFL=5.0" with linesp,\
     "../output/uFCT_Gal_BE_CFL10.txt" using 1:2 linecolor 2 pointtype 8\
       linetype 7 title "CFL=10.0" with linesp
