# note: this script requires Gnuplot version 4.6 or higher
# usage: gnuplot -e 'problem_ID=<problem ID>;
#           timeintegrator=<time integrator string>' solutions.gp

# list of possible input files to plot and their corresponding titles
file_initial  = "solution_".problem_ID."_initial"
file_exact    = "solution_".problem_ID."_exact"
file_galerkin = "solution_".problem_ID."_galerkin_".timeintegrator
file_low      = "solution_".problem_ID."_low_order_".timeintegrator
file_high     = "solution_".problem_ID."_high_order_".timeintegrator
file_EVFCT    = "solution_".problem_ID."_EVFCT_"     .timeintegrator
file_GalFCT   = "solution_".problem_ID."_GalFCT_"    .timeintegrator
file_DMPmin   = "DMPmin"
file_DMPmax   = "DMPmax"
file_list = file_initial." ".\
            file_exact." ".\
            file_galerkin." ".\
            file_low." ".\
            file_high." ".\
            file_EVFCT." ".\
            file_GalFCT." ".\
            file_DMPmin." ".\
            file_DMPmax
title_list = "Initial\
              Exact\
              Galerkin\
              Low-Order\
              High-Order\
              EV-FCT\
              Galerkin-FCT\
              DMP-min\
              DMP-max"
linetypes = "1 1 1 1 1 1 1 2 2"
linecolors = "7 -1 1 2 3 4 5 -1 -1"
symboltypes = "-1 -2 1 2 3 4 6 -2 -2"

# define is_missing(x) function for determining if an input file exists
outdir = "../output/problem_".problem_ID."/"
is_missing_aux(x)=system("ismissing.sh ".outdir.x)
is_missing(x)=int(is_missing_aux(x)+0)

# get a list of the possible input files that exist and their titles
existing_file_list  = ""
existing_title_list = ""
existing_lt_list = ""
existing_lc_list = ""
existing_sym_list = ""
do for [i=1:words(file_list)] {
   myfile = word(file_list,i).".gpl"
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
output_file = "../plots/solutions_".problem_ID."_".timeintegrator.".pdf"
set output '| ps2pdf - '.output_file
set ylabel "Angular Flux"
set xlabel "x"
set key top right

plot for [i=1:words(existing_file_list)] outdir.word(existing_file_list,i)\
   using 1:2 with linesp linetype word(existing_lt_list,i)\
   linecolor word(existing_lc_list,i)\
   pointtype word(existing_sym_list,i)\
   title word(existing_title_list,i)
