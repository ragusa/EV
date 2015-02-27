# note: this script requires Gnuplot version 4.6 or higher
# usage: gnuplot -e 'problem_ID=<problem ID>' solutions.gp

# list of possible input files to plot and their corresponding titles
file_list = "initial_solution\
             exact_solution\
             solution_galerkin\
             solution_low_order\
             solution_high_order\
             solution_FCT\
             min_values\
             max_values"
title_list = "Initial\
              Exact\
              Galerkin\
              Low-Order\
              High-Order\
              FCT\
              Min\
              Max"
linetypes = "1 1 1 1 1 1 2 2"
linecolors = "5 -1 1 3 4 2 -1 -1"
symboltypes = "-1 -2 1 2 3 4 -1 -1"


# define is_missing(x) function for determining if an input file exists
is_missing_aux(x)=system("ismissing.sh ".x)
is_missing(x)=int(is_missing_aux(x)+0)

# get a list of the possible input files that exist and their titles
existing_file_list  = ""
existing_title_list = ""
existing_lt_list = ""
existing_lc_list = ""
existing_sym_list = ""
do for [i=1:words(file_list)] {
   myfile = word(file_list,i)."_".problem_ID.".gpl"
   mytitle = word(title_list,i)
   mylt = word(linetypes,i)
   mylc = word(linecolors,i)
   mysym = word(symboltypes,i)
   if (!is_missing(myfile)) {
      existing_file_list = existing_file_list." ".myfile
      existing_title_list = existing_title_list." ".mytitle
      existing_lt_list = existing_lt_list." ".mylt
      existing_lc_list = existing_lc_list." ".mylc
      existing_sym_list = existing_sym_list." ".mysym
   }
}

set terminal postscript enhanced color
output_file = "../plots/solutions_".problem_ID.".pdf"
set output '| ps2pdf - '.output_file
set ylabel "Angular Flux"
set xlabel "x"
set key top right

plot for [i=1:words(existing_file_list)] "../output/".word(existing_file_list,i)\
   using 1:2 with linesp linetype word(existing_lt_list,i)\
   linecolor word(existing_lc_list,i)\
   pointtype word(existing_sym_list,i)\
   title word(existing_title_list,i)
