# submits a new job

import sys
import os
import subprocess
from subprocess import call

# check that correct number of command-line arguments was supplied
if len(sys.argv) != 3:
  sys.stderr.write("Usage : python %s <refinement> <n_processors>\n" % sys.argv[0])
  raise SystemExit(1)

# refinement level
r = int(sys.argv[1])
if r > 10:
  sys.stderr.write("Max refinement level is 10\n")
  raise SystemExit(1)
rstring = format(r,"02d")

# number of processors
p = int(sys.argv[2])
if p > 20:
  sys.stderr.write("Max number of cores is 20\n")
  raise SystemExit(1)
pstring = format(p,"02d")

# determine name of input file
input_file = "input_r" + rstring

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

# determine run number
runno = 0 # initialize run number to 0
job_not_found = True
newjob = ""
while job_not_found:
  runstring = format(runno,"02d")
  job = "p" + pstring + "_r" + rstring + "_" + runstring
  if job in joblist:
    if runno == 99:
      sys.stderr.write("There are already 100 runs for this study\n")
      raise SystemExit(1)
    runno = runno + 1
  else:
    newjob = job
    job_not_found = False

# form bsub command
outfile = "output/" + job + ".out"
#errfile = "output/" + job + ".err"
errfile = "output/err"
bsubcommand = str()
if p == 0: # p = 0 corresponds to serial
  bsubcommand = "bsub -c 100 -o " + outfile + " -e " + errfile + " -J " + newjob + \
                " transport_serial " + input_file
else:
  bsubcommand = "bsub -n " + str(p) + " -c 100 -o " + outfile + " -e " + errfile + " -J " + newjob + \
                " openmpi-mpirun transport_parallel " + input_file
os.system(bsubcommand)
