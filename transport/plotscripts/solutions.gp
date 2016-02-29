# note: this script requires Gnuplot version 4.6 or higher
# usage: gnuplot -e 'problem_name=<problem name>;
#           timeintegrator=<time integrator string>' solutions.gp

# list of possible input files to plot and their corresponding titles
file_initial  = "solution_initial"
file_exact    = "solution_exact"
file_galerkin = "solution_Gal_"   .timeintegrator
file_low      = "solution_low_"   .timeintegrator
file_high     = "solution_EV_"    .timeintegrator
file_EVFCT    = "solution_EVFCT_" .timeintegrator
file_GalFCT   = "solution_GalFCT_".timeintegrator
file_DMPmin  = "DMPmin"
file_DMPmax  = "DMPmax"
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
              EV\
              EV-FCT\
              Galerkin-FCT\
              DMP-min\
              DMP-max"
linetypes = "2 1 1 1 1 1 1 2 2"
linecolors = "7 -1 1 2 3 4 5 -1 -1"
symboltypes = "-1 -2 1 2 3 4 6 10 8"

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

      # check number of data points and do not use symbols if too many
      stats outdir.myfile using 1 noout
      n_data = STATS_records
      # print "Number of data points in ",myfile,": ",n_data
      #if (n_data > 200) {
      #   existing_sym_list   = existing_sym_list  ." -2"
      #} else {
         existing_sym_list   = existing_sym_list  ." ".mysym
      #}
   }
}

set terminal postscript enhanced color
output_file = "../output/".problem_name."/solutions_".timeintegrator.".pdf"
set output '| ps2pdf - '.output_file
set ylabel "Solution"
set xlabel "x"
set key top right

plot for [i=1:words(existing_file_list)] outdir.word(existing_file_list,i)\
   using 1:2 with linesp linetype word(existing_lt_list,i)\
   linecolor word(existing_lc_list,i)\
   pointtype word(existing_sym_list,i)\
   title word(existing_title_list,i)
