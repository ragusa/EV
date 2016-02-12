# create strong scaling plot data

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
  # open plot data output file
  output = open("strong_" + section + ".gpl",'w')
    
  # loop over number of processes
  mean_base = float()
  std_base = float()
  N_base = int()
  found_base = False
  for p in plist:
    case = "p" + format(p,"02d") + "_r" + format(r,"02d") + " " + section
    case_found = False
    for line in input:
      if case in line:
        words = line.split()
        mean = float(words[2])
        std = float(words[3])
        N = int(words[4])
        # use p = 1 as base
        if p == 1:
          mean_base = mean
          std_base = std
          N_base = N
          base_found = True
        # ensure that base was already found
        if not base_found:
          sys.stderr.write("Base has not yet been found\n")
          raise SystemExit(1)
        # check that the same number of runs were used for serial and parallel
        if N != N_base:
          sys.stderr.write("Different numbers of runs for case " + case + "\n")
          raise SystemExit(1)
        case_found = True
    if not case_found:
      sys.stderr.write("Data for case " + case + " was not found\n")
      raise SystemExit(1)
    # compute speedup and uncertainty bounds
    strong = mean_base / mean
    std = strong * sqrt((std/mean)**2 + (std_base/mean_base)**2)
    k = 1.96 * std / sqrt(N)
    # write line to output file
    output.write(str(p) + " " + str(strong) + " " + str(k) + "\n")

  # close plot data file
  output.close()
