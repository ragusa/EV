output_dir = "../output/mms_sinx_t/"
output_file = output_dir."solutions.pdf"

# create plot
set terminal postscript enhanced color
set output '| ps2pdf - '.output_file
set ylabel "Solution"
set xlabel "x"
file_list = "solution_exact.gpl solution_EV_SSP3.gpl"
symbol_list = "2 3"
#plot for [i=1:words(file_list)] output_dir.word(file_list,i)\
#  u 1:2 w linesp lt 1\
#  pt word(symbol_list,i)+"0"
plot output_dir."solution_exact.gpl" using 1:2 title "Exact" with linesp pointtype int("2"),\
     output_dir."solution_EV_SSP3.gpl" using 1:2 title "EV" with linesp pointtype int("3")
