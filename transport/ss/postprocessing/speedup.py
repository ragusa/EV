# create speedup plot data

import sys
from math import sqrt

# refinement level used for speedup study
r = 9
# list of number of processes
plist = list([1,2,4,8,16,20])
# list of timed sections
sectionlist = list(["assemble","initialize","output","setup","solve"])

# get lines from input file
inputfile = open("results.out",'r')
input = inputfile.readlines()
inputfile.close()

# loop over sections
for section in sectionlist:
  # get serial data
  case = "p00_r" + format(r,"02d") + " " + section
  mean_serial = float()
  std_serial = float()
  N_serial = int()
  case_found = False
  for line in input:
    if case in line:
      words = line.split()
      mean_serial = float(words[2])
      std_serial = float(words[3])
      N_serial = int(words[4])
      case_found = True
  if not case_found:
    sys.stderr.write("Data for case " + case + " was not found\n")
    raise SystemExit(1)

  # open plot data output file
  output = open("speedup_" + section + ".gpl",'w')
    
  # loop over number of processes
  for p in plist:
    case = "p" + format(p,"02d") + "_r" + format(r,"02d") + " " + section
    case_found = False
    for line in input:
      if case in line:
        words = line.split()
        mean = float(words[2])
        std = float(words[3])
        N = int(words[4])
        # check that the same number of runs were used for serial and parallel
        if N != N_serial:
          sys.stderr.write("Different numbers of runs for case " + case + "\n")
          raise SystemExit(1)
        case_found = True
    if not case_found:
      sys.stderr.write("Data for case " + case + " was not found\n")
      raise SystemExit(1)
    # compute speedup and uncertainty bounds
    speedup = mean_serial / mean
    std = speedup * sqrt((std/mean)**2 + (std_serial/mean_serial)**2)
    k = 1.96 * std / sqrt(N)
    # write line to output file
    output.write(str(p) + " " + str(speedup) + " " + str(k) + "\n")

  # close plot data file
  output.close()
