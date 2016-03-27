# note: this script requires Gnuplot version 4.6 or higher
# usage: gnuplot -e 'parameter_name=<parameter name>;
#   parameter_values=<parameter values>' study_solution.gp

filebase = "solution"

# list of possible input files to plot and their corresponding titles
file_exact = filebase."_exact"
file_0     = filebase."0"
file_1     = filebase."1"
file_2     = filebase."2"
file_3     = filebase."3"
file_4     = filebase."4"
file_5     = filebase."5"
file_6     = filebase."6"
file_7     = filebase."7"
file_8     = filebase."8"
file_9     = filebase."9"
file_list = file_exact." ".\
            file_0." ".\
            file_1." ".\
            file_2." ".\
            file_3." ".\
            file_4." ".\
            file_5." ".\
            file_6." ".\
            file_7." ".\
            file_8." ".\
            file_9
linetypes = "1 1 1 1 1 1 1 1 1 1 1"
linecolors = "-1 0 1 2 3 4 5 6 7 8 9"
symboltypes = "-2 0 1 2 3 4 5 6 7 8 9"

# define is_missing(x) function for determining if an input file exists
outdir = "./output/"
is_missing_aux(x)=system("../../plotscripts/ismissing.sh ".outdir.x)
is_missing(x)=int(is_missing_aux(x)+0)

# set print file to STDOUT
set print "-"

# get a list of the possible input files that exist and their titles
existing_file_list  = ""
existing_title_list = ""
existing_lt_list = ""
existing_lc_list = ""
existing_sym_list = ""
do for [i=1:words(file_list)] {
   myfile  = word(file_list,i).".gpl"
   mytitle = "blah"
   if (i == 1) {
     mytitle = "Exact"
   } else {
     mytitle = parameter_name."=".word(parameter_values,i-1)
   }
   mylt    = word(linetypes,i)
   mylc    = word(linecolors,i)
   mysym   = word(symboltypes,i)
   if (!is_missing(myfile)) {
      existing_file_list  = existing_file_list ." ".myfile
      existing_title_list = existing_title_list." ".mytitle
      existing_lt_list    = existing_lt_list   ." ".mylt
      existing_lc_list    = existing_lc_list   ." ".mylc

      # check number of data points and do not use symbols if too many
      stats outdir.myfile using 1 noout
      n_data = STATS_records
      # print "Number of data points in ",myfile,": ",n_data
      if (n_data > 200) {
         existing_sym_list   = existing_sym_list  ." -2"
      } else {
         existing_sym_list   = existing_sym_list  ." ".mysym
      }
   }
}

set terminal postscript enhanced color
output_file = "solutions.pdf"
set output '| ps2pdf - '.output_file
set ylabel "Solution"
set xlabel "x"
set key top right

plot for [i=1:words(existing_file_list)] outdir.word(existing_file_list,i)\
  using 1:2 with linesp linetype word(existing_lt_list,i)\
  linecolor word(existing_lc_list,i)\
  pointtype word(existing_sym_list,i)\
  title word(existing_title_list,i)
