# note: this script requires Gnuplot version 4.6 or higher
# usage: gnuplot -e 'problem_ID=<problem ID>;
#           timeintegrator=<time integrator string>' solutions.gp

# list of possible input files to plot and their corresponding titles
file_initial  = "solution_".problem_ID."_initial"
file_exact    = "solution_".problem_ID."_exact"
file_galerkin = "solution_".problem_ID."_Gal_".timeintegrator
file_low      = "solution_".problem_ID."_low_".timeintegrator
file_high     = "solution_".problem_ID."_EV_".timeintegrator
file_EVFCT    = "solution_".problem_ID."_EVFCT_"     .timeintegrator
file_GalFCT   = "solution_".problem_ID."_GalFCT_"    .timeintegrator
file_LowDMPmin = "DMPmin_Low"
file_LowDMPmax = "DMPmax_Low"
file_GalFCTDMPmin = "DMPmin_GalFCT"
file_GalFCTDMPmax = "DMPmax_GalFCT"
file_EVFCTDMPmin  = "DMPmin_EVFCT"
file_EVFCTDMPmax  = "DMPmax_EVFCT"
file_list = file_initial." ".\
            file_exact." ".\
            file_galerkin." ".\
            file_low." ".\
            file_high." ".\
            file_EVFCT." ".\
            file_GalFCT." ".\
            file_LowDMPmin." ".\
            file_LowDMPmax." ".\
            file_GalFCTDMPmin." ".\
            file_GalFCTDMPmax." ".\
            file_EVFCTDMPmin." ".\
            file_EVFCTDMPmax
title_list = "Initial\
              Exact\
              Galerkin\
              Low-Order\
              EV\
              EV-FCT\
              Galerkin-FCT\
              DMP-min\
              DMP-max\
              DMP-min-Gal-FCT\
              DMP-max-Gal-FCT\
              DMP-min-EV-FCT\
              DMP-max-EV-FCT"
linetypes = "1 1 1 1 1 1 1 2 2 2 2 4 4"
linecolors = "7 -1 1 2 3 4 5 -1 -1 -1 -1 -1 -1"
symboltypes = "-1 -2 1 2 3 4 6 -2 -2 -2 -2 -2 -2"

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
set ylabel "Solution"
set xlabel "x"
set key top right

plot for [i=1:words(existing_file_list)] outdir.word(existing_file_list,i)\
   using 1:2 with linesp linetype word(existing_lt_list,i)\
   linecolor word(existing_lc_list,i)\
   pointtype word(existing_sym_list,i)\
   title word(existing_title_list,i)
