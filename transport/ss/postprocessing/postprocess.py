# postprocess raw output files

import numpy as np
import os
import sys
from math import sqrt

# list of timed sections
sectionlist = list(["assemble","initialize","output","setup","solve"])

# list of files in output directory
filelist = os.listdir("../output")
if "err" in filelist:
  filelist.remove("err")

# create dictionary with runcase -> casefilelist
casefiledict = dict()
for file in filelist:
  # extract case name from file name
  fields = (file.split(".out")[0]).split("_")
  case = fields[0] + "_" + fields[1]
  if case in casefiledict:
    # add file to existing file list
    casefiledict[case].append(file)
  else:
    # create new file list for case
    casefiledict[case] = list([file])

# open results file to write
fresults = open("results.out",'w')

# for each case, loop through its files and extract times
for case in casefiledict.keys():
  # create dictionary for list of times for each section
  times = dict()
  for section in sectionlist:
    times[section] = list()
  # loop over files for the case
  for file in casefiledict[case]:
    # create dictionary for flags for each section being found
    sectionfound = dict()
    for section in sectionlist:
      sectionfound[section] = False
    # open file
    f = open("../output/" + file,'r')
    # loop over lines of file
    for line in f:
      # loop over sections
      for section in sectionlist:
        # check if line contains time for section
        if "| " + section in line:
          # get time
          entries = line.split("|")
          entry = entries[3].lstrip()
          timestring = entry.split("s")[0]
          time = float(timestring)
          # add time to section time list
          times[section].append(time)
          # mark section as being found
          sectionfound[section] = True
    # close file
    f.close()
    # check that all sections were found in file
    for section in sectionlist:
      if not sectionfound[section]:
        sys.stderr.write("Section \"" + section + "\" not found in " + file + "\n")
        raise SystemExit(1)
  # calculate statistics on times for case
  for section in sectionlist:
    # number of experiments
    n = len(times[section])
    # compute mean
    mean = np.mean(times[section])
    # compute standard deviation
    std = np.std(times[section],ddof=1)
    # write to file
    fresults.write(case + " " + section + " " + str(mean) + " " +
      str(std) + " " + str(n) + "\n")

# close results file
fresults.close()
