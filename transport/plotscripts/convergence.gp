# note: this script requires Gnuplot version 4.6 or higher
# usage: gnuplot -e 'problem_ID=<problem ID>;
#           timeintegrator=<time integrator string>' convergence.gp

# Reference slopes
slope1 = 0.25
slope2 = 0.50
slope3 = 0.75

# Convergence mode; 1 = space, 2 = time
conv_mode = 1

# Steady-state mode; 1 = transient, 2 = steady-state
ss_mode = 1

# list of possible input files to plot and their corresponding titles
file_galerkin = "convergence_".problem_ID."_galerkin_".timeintegrator
file_low      = "convergence_".problem_ID."_low_order_".timeintegrator
file_high     = "convergence_".problem_ID."_high_order_".timeintegrator
file_EVFCT    = "convergence_".problem_ID."_EVFCT_"     .timeintegrator
file_GalFCT   = "convergence_".problem_ID."_GalFCT_"    .timeintegrator
file_list = file_galerkin." ".\
            file_low." ".\
            file_high." ".\
            file_EVFCT." ".\
            file_GalFCT
title_list = "Galerkin\
              Low-Order\
              High-Order\
              EV-FCT\
              Galerkin-FCT"
linetypes = "1 1 1 1 1"
linecolors = "1 2 3 4 5"
symboltypes = "1 2 3 4 6"

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

# determine h_min based on convergence mode and then get minimum dx or dt
if (conv_mode == 1) \
  h_col = 4; \
else \
  h_col = 5
stats outdir.word(existing_file_list,words(existing_file_list)) using h_col noout
h_min = STATS_min

# get minimum L1 error
if (ss_mode == 1) \
  L1_col = 7; \
else \
  L1_col = 5
stats outdir.word(existing_file_list,words(existing_file_list)) using L1_col noout
L1_min = STATS_min
L1_max = STATS_max

# get minimum L2 error
if (ss_mode == 1) \
  L2_col = 9; \
else \
  L2_col = 7
stats outdir.word(existing_file_list,words(existing_file_list)) using L2_col noout
L2_min = STATS_min
L2_max = STATS_max

# define reference slope functions
L1c1 = (1.0/h_min)**slope1 * L1_min
L1c2 = (1.0/h_min)**slope2 * L1_min
L1c3 = (1.0/h_min)**slope3 * L1_min
L1ref1(x) = L1c1 * x**slope1
L1ref2(x) = L1c2 * x**slope2
L1ref3(x) = L1c3 * x**slope3
L1string1 = sprintf("m=%.2f slope",slope1)
L1string2 = sprintf("m=%.2f slope",slope2)
L1string3 = sprintf("m=%.2f slope",slope3)

L2c1 = (1.0/h_min)**slope1 * L2_min
L2c2 = (1.0/h_min)**slope2 * L2_min
L2c3 = (1.0/h_min)**slope3 * L2_min
L2ref1(x) = L2c1 * x**slope1
L2ref2(x) = L2c2 * x**slope2
L2ref3(x) = L2c3 * x**slope3
L2string1 = sprintf("m=%.2f slope",slope1)
L2string2 = sprintf("m=%.2f slope",slope2)
L2string3 = sprintf("m=%.2f slope",slope3)

set terminal postscript enhanced color
if (conv_mode == 1) \
  set xlabel "Mesh Size"; \
else \
  set xlabel "Time Step Size"
set logscale xy
set key top left
set format y "10^{%L}"
set format x "10^{%L}"

output_file = "../plots/convergence_".problem_ID."_".timeintegrator.".pdf"
set output '| ps2pdf - '.output_file

# create multiplot
set multiplot layout 1, 2

# Plot L-1 error
if (L1_max > 1) \
  set yrange[*:1]; \
else \
  set yrange[*:*]
set ylabel "L-1 Error"
plot for [i=1:words(existing_file_list)] outdir.word(existing_file_list,i)\
   using h_col:L1_col with linesp linetype word(existing_lt_list,i)\
   linecolor word(existing_lc_list,i)\
   pointtype word(existing_sym_list,i)\
   title word(existing_title_list,i),\
   L1ref1(x) title L1string1 with lines linestyle 2 linecolor 9,\
   L1ref2(x) title L1string2 with lines linestyle 2 linecolor 9,\
   L1ref3(x) title L1string3 with lines linestyle 2 linecolor 9

# Plot L-2 error
if (L2_max > 1) \
  set yrange[*:1]; \
else \
  set yrange[*:*]
set ylabel "L-2 Error"
plot for [i=1:words(existing_file_list)] outdir.word(existing_file_list,i)\
   using h_col:L2_col with linesp linetype word(existing_lt_list,i)\
   linecolor word(existing_lc_list,i)\
   pointtype word(existing_sym_list,i)\
   title word(existing_title_list,i),\
   L2ref1(x) title L2string1 with lines linestyle 2 linecolor 9,\
   L2ref2(x) title L2string2 with lines linestyle 2 linecolor 9,\
   L2ref3(x) title L2string3 with lines linestyle 2 linecolor 9

unset multiplot
