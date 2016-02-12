set terminal postscript enhanced color
output_file = "plot2d.pdf"
set output output_file
input_file = "limiter.txt"

set xlabel "i"
set ylabel "j"
n_rows = `awk 'END {print NR}' limiter.txt`
set xrange[0:n_rows-1]
set yrange[0:n_rows-1]
set key outside
plot input_file matrix with image notitle
