# print summary of currently running and already run jobs

import sys
import os
import subprocess
from subprocess import call

# create joblist
joblist = list()

# get list of jobs in queue
subp = subprocess.Popen("bjobs -w", shell=True, stdout=subprocess.PIPE)
subp_out = subp.stdout.read()
lines = subp_out.split("\n")
for line in lines:
  if line: # exclude last line, which is empty
    words = line.split()
    job = words[6]
    joblist.append(job)

# get list of jobs already run
filelist = os.listdir("output")
for file in filelist:
  dotsep = file.split(".")
  job = dotsep[0]
  joblist.append(job)

# loop over jobs and add to counts
casedict = dict()
for job in joblist:
  if job != "err" and job != "JOB_NAME":
    fields = job.split("_")
    case = fields[0] + "_" + fields[1]
    if case in casedict:
      casedict[case] = casedict[case] + 1
    else:
      casedict[case] = 1

# get sorted list of keys
sortedkeys = sorted(casedict.keys())

# print header
print("p  r  runs");
print("----------");

# loop over cases
for key in sortedkeys:
  keyfields = key.split("_")
  pstring = keyfields[0].split("p")[1]
  rstring = keyfields[1].split("r")[1]
  print(pstring + " " + rstring + " " + str(casedict[key]))

