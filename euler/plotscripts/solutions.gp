# note: this script requires Gnuplot version 4.6 or higher
# usage: gnuplot -e 'problem_name=<problem name>;
#           timeintegrator=<time integrator string>;
#           component=<component of solution>;
#           component_filename=<filename for component>;
#           component_column=<column of component>' solutions.gp

# list of possible input files to plot and their corresponding titles
file_initial  = "solution_initial"
file_exact    = "solution_exact"
file_galerkin = "solution_Gal_"   .timeintegrator
file_low      = "solution_low_"   .timeintegrator
file_high     = "solution_EV_"    .timeintegrator
file_EVFCT    = "solution_EVFCT_" .timeintegrator
file_GalFCT   = "solution_GalFCT_".timeintegrator
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
outdir = "../output/".problem_name."/"
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
output_file = outdir.component_filename."_".timeintegrator.".pdf"
set output '| ps2pdf - '.output_file
set ylabel component
set xlabel "x"
set key top right

plot for [i=1:words(existing_file_list)] outdir.word(existing_file_list,i)\
   using 1:int(component_column) with linesp linetype word(existing_lt_list,i)\
   linecolor word(existing_lc_list,i)\
   pointtype word(existing_sym_list,i)\
   title word(existing_title_list,i)
