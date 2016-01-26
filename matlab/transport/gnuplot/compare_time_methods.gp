output_file = "../output/compare_time_methods.pdf"

# create plot
set terminal pdfcairo
set output output_file
set ylabel "Solution"
set xlabel "x"
plot "../output/uexact.txt" using 1:2 linecolor 0 pointtype 0 title "Exact" with linesp,\
     "../output/uFCT_Gal_BE.txt" using 1:2 linecolor 1 pointtype 1 title "BE" with linesp,\
     "../output/uFCT_Gal_FE.txt" using 1:2 linecolor 2 pointtype 2 title "FE" with linesp,\
     "../output/uFCT_Gal_SSPRK33.txt" using 1:2 linecolor 3 pointtype 3 title "SSPRK33" with linesp
