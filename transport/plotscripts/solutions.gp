# note: this script requires Gnuplot version 4.6 or higher
# usage: gnuplot -e 'problem_ID=<problem ID>' solutions.gp

# list of possible input files to plot and their corresponding titles
file_list = "exact_solution\
             solution_none\
             solution_old_first_order\
             solution_old_entropy\
             solution_low_order\
             solution_high_order"
title_list = "Exact\
              No-Viscosity\
              Old-First-Order\
              Old-Entropy\
              Low-Order\
              High-Order"

# define is_missing(x) function for determining if an input file exists
is_missing_aux(x)=system("ismissing.sh ".x)
is_missing(x)=int(is_missing_aux(x)+0)

# get a list of the possible input files that exist and their titles
existing_file_list  = ""
existing_title_list = ""
do for [i=1:words(file_list)] {
   myfile = word(file_list,i)."_".problem_ID.".gpl"
   mytitle = word(title_list,i)
   if (!is_missing(myfile)) {
      existing_file_list = existing_file_list." ".myfile
      existing_title_list = existing_title_list." ".mytitle
   }
}

set terminal postscript enhanced color
output_file = "solutions_".problem_ID.".pdf"
set output '| ps2pdf - '.output_file
set ylabel "Angular Flux"
set xlabel "x"
set key bottom right

plot for [i=1:words(existing_file_list)] "../output/".word(existing_file_list,i)\
   using 1:2 with lines title word(existing_title_list,i)
