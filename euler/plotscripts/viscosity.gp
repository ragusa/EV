# note: this script requires Gnuplot version 4.6 or higher
# usage: gnuplot -e 'problem_name=<problem name>;
#           timeintegrator=<time integrator string>'
#           viscosity.gp

filebase = "viscosity"

# list of possible input files to plot and their corresponding titles
file_loworder  = "loworderviscosity_" .timeintegrator
file_entropy   = "entropyviscosity_"  .timeintegrator
file_highorder = "highorderviscosity_".timeintegrator
file_list = file_loworder." ".\
            file_entropy." ".\
            file_highorder." "
title_list = "Low-Order\
              Entropy\
              High-Order"
linetypes = "1 1 1"
linecolors = "6 3 2"
symboltypes = "-2 -2 -2"

# define is_missing(x) function for determining if an input file exists
outdir = "../output/".problem_name."/"
is_missing_aux(x)=system("ismissing.sh ".outdir.x)
is_missing(x)=int(is_missing_aux(x)+0)

# set print file to STDOUT
set print "-"

# get a list of the possible input files that exist and their titles
existing_file_list  = ""
existing_title_list = ""
existing_lt_list = ""
existing_lc_list = ""
existing_sym_list = ""
do for [i=1:words(file_list)] {
   myfile  = word(file_list,i).".gpl"
   mytitle = word(title_list,i)
   mylt    = word(linetypes,i)
   mylc    = word(linecolors,i)
   mysym   = word(symboltypes,i)
   if (!is_missing(myfile)) {
      existing_file_list  = existing_file_list ." ".myfile
      existing_title_list = existing_title_list." ".mytitle
      existing_lt_list    = existing_lt_list   ." ".mylt
      existing_lc_list    = existing_lc_list   ." ".mylc
      existing_sym_list   = existing_sym_list  ." ".mysym
   }
}

set terminal pdfcairo
output_file = outdir."viscosity_".timeintegrator.".pdf"
set output output_file
set ylabel "Viscosity"
set xlabel "x"
#set logscale y
set key top right
set style fill transparent solid 0.5 noborder

plot for [i=1:words(existing_file_list)]\
  outdir.word(existing_file_list,i) using 1:2\
  with filledcurves y1=0\
  linetype word(existing_lt_list,i)\
  linecolor word(existing_lc_list,i)\
  title word(existing_title_list,i)

