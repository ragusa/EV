# define is_missing(x) function to determine if input file exists
is_missing_aux(x)=system("ismissing.sh ".x)
is_missing(x)=int(is_missing_aux(x)+0)

set terminal postscript enhanced color
set output '| ps2pdf - test_problem_1.pdf'
set ylabel "Angular Flux"
set xlabel "x"
set key bottom right

file_list = "analytic_solution_test_problem_1.dat\
             solution_none_1.gpl\
             solution_first_order_1.gpl\
             solution_entropy_1.gpl\
             solution_max_principle_1.gpl"
title_list = "Analytic\
              No-Viscosity\
              Old-First-Order\
              Old-Entropy\
              First-Order"

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

plot for [i=1:words(existing_file_list)] "../output/".word(existing_file_list,i)\
   using 1:2 with lines title word(existing_title_list,i)

reset
set output
set terminal pop
