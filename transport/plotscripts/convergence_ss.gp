# note: this script requires Gnuplot version 4.6 or higher
# usage: gnuplot -e 'problem_ID=<problem ID>' convergence.gp

# list of possible input files to plot and their corresponding titles
file_list = "convergence_galerkin\
             convergence_low_order\
             convergence_high_order\
             convergence_FCT"
title_list = "Galerkin\
              Low-Order\
              High-Order\
              FCT"
linetypes = "1 1 1 1"
linecolors = "1 3 4 2"
symboltypes = "1 2 3 4"


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

stats "../output/".word(existing_file_list,words(existing_file_list)) using 4 noout
h_min = STATS_min

stats "../output/".word(existing_file_list,words(existing_file_list)) using 5 noout
L1_min = STATS_min

stats "../output/".word(existing_file_list,words(existing_file_list)) using 7 noout
L2_min = STATS_min

# define reference slope functions
L1c1 = (1.0/h_min)**1 * L1_min
L1c2 = (1.0/h_min)**2 * L1_min
L1c3 = (1.0/h_min)**3 * L1_min
L1ref1(x) = L1c1 * x**1   
L1ref2(x) = L1c2 * x**2  
L1ref3(x) = L1c3 * x**3 

L2c1 = (1.0/h_min)**1 * L2_min
L2c2 = (1.0/h_min)**2 * L2_min
L2c3 = (1.0/h_min)**3 * L2_min
L2ref1(x) = L2c1 * x**1   
L2ref2(x) = L2c2 * x**2  
L2ref3(x) = L2c3 * x**3 

set terminal postscript enhanced color
set xlabel "Mesh Size"
set logscale xy
set key top right
set format y "10^{%L}"
set format x "10^{%L}"

output_file = "../plots/convergence_".problem_ID.".pdf"
set output '| ps2pdf - '.output_file

set multiplot layout 1, 2

# Plot L-1 error
#output_file = "../plots/convergence_".problem_ID."_L1.pdf"
#set output '| ps2pdf - '.output_file
set ylabel "L-1 Error"
plot for [i=1:words(existing_file_list)] "../output/".word(existing_file_list,i)\
   using 4:5 with linesp linetype word(existing_lt_list,i)\
   linecolor word(existing_lc_list,i)\
   pointtype word(existing_sym_list,i)\
   title word(existing_title_list,i),\
   L1ref1(x) title "m=1 slope" with lines linestyle 2 linecolor 9,\
   L1ref2(x) title "m=2 slope" with lines linestyle 2 linecolor 9,\
   L1ref3(x) title "m=3 slope" with lines linestyle 2 linecolor 9

# Plot L-2 error
#output_file = "../plots/convergence_".problem_ID."_L2.pdf"
#set output '| ps2pdf - '.output_file
set ylabel "L-2 Error"
plot for [i=1:words(existing_file_list)] "../output/".word(existing_file_list,i)\
   using 4:7 with linesp linetype word(existing_lt_list,i)\
   linecolor word(existing_lc_list,i)\
   pointtype word(existing_sym_list,i)\
   title word(existing_title_list,i),\
   L2ref1(x) title "m=1 slope" with lines linestyle 2 linecolor 9,\
   L2ref2(x) title "m=2 slope" with lines linestyle 2 linecolor 9,\
   L2ref3(x) title "m=3 slope" with lines linestyle 2 linecolor 9

unset multiplot
