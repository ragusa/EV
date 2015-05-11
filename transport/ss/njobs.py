# count number of jobs in queue

import subprocess

subp = subprocess.Popen("bjobs -w", shell=True, stdout=subprocess.PIPE)
subp_out = subp.stdout.read()
lines = subp_out.split("\n")
njobs = 0
for line in lines:
  if line: # exclude last line, which is empty
    if line.split()[0] != "JOBID":
      njobs = njobs + 1
print("Number of jobs = " + str(njobs))
