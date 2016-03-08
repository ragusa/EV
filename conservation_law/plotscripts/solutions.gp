# note: this script requires Gnuplot version 4.6 or higher
# usage: gnuplot -e 'problem_name=<problem name>;
#           timeintegrator=<time integrator string>;
#           quantity=<quantity>;
#           column=<column of quantity in input file>' solutions.gp

filebase = "solution"

# list of possible input files to plot and their corresponding titles
file_initial  = filebase."_initial"
file_exact    = filebase."_exact"
file_Gal      = filebase."_Gal_"   .timeintegrator
file_Lax      = filebase."_Lax_"   .timeintegrator
file_DMP      = filebase."_DMP_"   .timeintegrator
file_DIV      = filebase."_DIV_"   .timeintegrator
file_DID      = filebase."_DID_"   .timeintegrator
file_EV       = filebase."_EV_"    .timeintegrator
file_EVFCT    = filebase."_EVFCT_" .timeintegrator
file_GalFCT   = filebase."_GalFCT_".timeintegrator
file_DMPmin = "DMPmin"
file_DMPmax = "DMPmax"
file_GalFCTmin = "GalFCTmin"
file_GalFCTmax = "GalFCTmax"
file_EVFCTmin  = "EVFCTmin"
file_EVFCTmax  = "EVFCTmax"
file_list = file_initial." ".\
            file_exact." ".\
            file_Gal." ".\
            file_Lax." ".\
            file_DMP." ".\
            file_DIV." ".\
            file_DID." ".\
            file_EV." ".\
            file_EVFCT." ".\
            file_GalFCT." ".\
            file_DMPmin." ".\
            file_DMPmax." ".\
            file_GalFCTmin." ".\
            file_GalFCTmax." ".\
            file_EVFCTmin." ".\
            file_EVFCTmax
title_list = "Initial\
              Exact\
              Galerkin\
              Lax\
              DMP\
              DI-Viscosity\
              DI-Diffusion\
              EV\
              EV-FCT\
              Galerkin-FCT\
              DMP-min\
              DMP-max\
              DMP-min-Gal-FCT\
              DMP-max-Gal-FCT\
              DMP-min-EV-FCT\
              DMP-max-EV-FCT"
linetypes = "2 1 1 2 4 1 3 1 1 1 2 2 2 2 4 4"
linecolors = "-1 -1 1 2 2 2 2 3 4 5 -1 -1 -1 -1 -1 -1"
symboltypes = "-2 -2 1 4 3 2 1 3 4 6 -2 -2 -2 -2 -2 -2"

# define is_missing(x) function for determining if an input file exists
outdir = "../output/".problem_name."/"
is_missing_aux(x)=system("ismissing.sh ".outdir.x)
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
   mytitle = word(title_list,i)
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

# determine y label
if (quantity eq "velocity") {
   quantity_ylabel = "Velocity"
} else { if (quantity eq "density") {
   quantity_ylabel = "Density"
} else { if (quantity eq "momentum") {
   quantity_ylabel = "Momentum"
} else { if (quantity eq "energy") {
   quantity_ylabel = "Energy"
} else { if (quantity eq "height") {
   quantity_ylabel = "Height"
} else { if (quantity eq "waterlevel") {
   quantity_ylabel = "Water Level"
} else { if (quantity eq "angularflux") {
   quantity_ylabel = "Angular Flux"
} else {
   quantity_ylabel = "Unknown"
}}}}}}}

set terminal postscript enhanced color
output_file = outdir.quantity."_".timeintegrator.".pdf"
set output '| ps2pdf - '.output_file
set ylabel quantity_ylabel
set xlabel "x"
set key top right

# if water level plot, then plot with bathymetry function
if (quantity eq "waterlevel") {

   bathymetry_file = file_galerkin.".gpl"
   if (is_missing(bathymetry_file)) {
      bathymetry_file = file_low.".gpl"
      if (is_missing(bathymetry_file)) {
         bathymetry_file = file_high.".gpl"
         if (is_missing(bathymetry_file)) {
            bathymetry_file = file_EVFCT.".gpl"
            if (is_missing(bathymetry_file)) {
               bathymetry_file = file_GalFCT.".gpl"
            }
         }
      }
   }

   set style fill pattern 7
   plot for [i=1:words(existing_file_list)] outdir.word(existing_file_list,i)\
      using 1:int(column) with linesp linetype word(existing_lt_list,i)\
      linecolor word(existing_lc_list,i)\
      pointtype word(existing_sym_list,i)\
      title word(existing_title_list,i),\
      outdir.bathymetry_file using 1:4 with filledcurves y1=0 linecolor 0\
        linetype 1 title "Bottom topography"

} else {
   plot for [i=1:words(existing_file_list)] outdir.word(existing_file_list,i)\
      using 1:int(column) with linesp linetype word(existing_lt_list,i)\
      linecolor word(existing_lc_list,i)\
      pointtype word(existing_sym_list,i)\
      title word(existing_title_list,i)
}
