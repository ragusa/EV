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

set terminal postscript enhanced color
set xlabel "Mesh Size"
set logscale xy
set key top right
set format y "10^{%L}"
set format x "10^{%L}"
#set yrange [:1e-1]

# define reference slope functions
L1c1 = 512**1 * 1.0e-4
L1c2 = 512**2 * 1.0e-4
L1c3 = 512**3 * 1.0e-4
L1ref1(x) = L1c1 * x**1   
L1ref2(x) = L1c2 * x**2  
L1ref3(x) = L1c3 * x**3 

L2c1 = 512**1 * 2.0e-4
L2c2 = 512**2 * 2.0e-4
L2c3 = 512**3 * 2.0e-4
L2ref1(x) = L2c1 * x**1   
L2ref2(x) = L2c2 * x**2  
L2ref3(x) = L2c3 * x**3 

# Plot L-1 error
output_file = "convergence_".problem_ID."_L1.pdf"
set output '| ps2pdf - '.output_file
set ylabel "L-1 Error"
plot for [i=1:words(existing_file_list)] "../output/".word(existing_file_list,i)\
   using 4:7 with linesp linetype word(existing_lt_list,i)\
   linecolor word(existing_lc_list,i)\
   pointtype word(existing_sym_list,i)\
   title word(existing_title_list,i),\
   L1ref1(x) title "m=1 slope" with lines linestyle 2 linecolor 9,\
   L1ref2(x) title "m=2 slope" with lines linestyle 2 linecolor 9,\
   L1ref3(x) title "m=3 slope" with lines linestyle 2 linecolor 9

# Plot L-2 error
output_file = "convergence_".problem_ID."_L2.pdf"
set output '| ps2pdf - '.output_file
set ylabel "L-2 Error"
plot for [i=1:words(existing_file_list)] "../output/".word(existing_file_list,i)\
   using 4:9 with linesp linetype word(existing_lt_list,i)\
   linecolor word(existing_lc_list,i)\
   pointtype word(existing_sym_list,i)\
   title word(existing_title_list,i),\
   L2ref1(x) title "m=1 slope" with lines linestyle 2 linecolor 9,\
   L2ref2(x) title "m=2 slope" with lines linestyle 2 linecolor 9,\
   L2ref3(x) title "m=3 slope" with lines linestyle 2 linecolor 9
