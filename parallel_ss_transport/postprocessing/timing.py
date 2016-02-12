# create timing plot data

import sys
from math import sqrt

# refinement levels used in plot
rlist = list([5,7,9])
# list of number of processes
plist = list([0,1,8,20])
# list of timed sections
sectionlist = list(["assemble","initialize","output","setup","solve"])

# get lines from input file
inputfile = open("results.out",'r')
input = inputfile.readlines()
inputfile.close()

# loop over sections
for section in sectionlist:
  # loop over number of processes
  for p in plist:
    # open plot data output file
    output = open("timing_p" + format(p,"02d") + "_" + section + ".gpl",'w')
    
    # loop over refinement levels
    for r in rlist:
      case = "p" + format(p,"02d") + "_r" + format(r,"02d") + " " + section
      case_found = False
      for line in input:
        if case in line:
          words = line.split()
          mean = float(words[2])
          std = float(words[3])
          N = int(words[4])
          case_found = True
      if not case_found:
        sys.stderr.write("Data for case " + case + " was not found\n")
        raise SystemExit(1)
      # compute speedup and uncertainty bounds
      k = 1.96 * std / sqrt(N)
      # compute size corresponding to refinement level
      size = (2**r + 1)**2
      # write line to output file
      output.write(str(size) + " " + str(mean) + " " + str(k) + "\n")

    # close plot data file
    output.close()
