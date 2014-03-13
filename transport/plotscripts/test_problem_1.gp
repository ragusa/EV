# note: this script requires Gnuplot version 4.6 or higher

# list of possible input files to plot and their corresponding titles
file_list = "exact_solution_1.gpl\
             solution_none_1.gpl\
             solution_old_first_order_1.gpl\
             solution_old_entropy_1.gpl\
             solution_first_order_1.gpl\
             solution_entropy_1.gpl"
title_list = "Exact\
              No-Viscosity\
              Old-First-Order\
              Old-Entropy\
              First-Order\
              Entropy"

# define is_missing(x) function for determining if an input file exists
is_missing_aux(x)=system("ismissing.sh ".x)
is_missing(x)=int(is_missing_aux(x)+0)

# get a list of the possible input files that exist and their titles
existing_file_list  = ""
existing_title_list = ""
do for [i=1:words(file_list)] {
   myfile = word(file_list,i)
   mytitle = word(title_list,i)
   if (!is_missing(myfile)) {
      existing_file_list = existing_file_list." ".myfile
      existing_title_list = existing_title_list." ".mytitle
   }
}

set terminal postscript enhanced color
set output '| ps2pdf - test_problem_1.pdf'
set ylabel "Angular Flux"
set xlabel "x"
set key bottom right

plot for [i=1:words(existing_file_list)] "../output/".word(existing_file_list,i)\
   using 1:2 with lines title word(existing_title_list,i)
