# create plot
set terminal pdfcairo
output_file = "limiter.pdf"
set output output_file
plot "limiter.txt" matrix with image
