# run multiple jobs

import sys
import os

# check that correct number of command-line arguments was supplied
if len(sys.argv) != 4:
  sys.stderr.write("Usage : python %s <refinement> <n_processors> <n_runs>\n" % sys.argv[0])
  raise SystemExit(1)

# refinement level
rstring = sys.argv[1]
r = int(rstring)
if r > 10:
  sys.stderr.write("Max refinement level is 10\n")
  raise SystemExit(1)

# number of processors
pstring = sys.argv[2]
p = int(pstring)
if p > 20:
  sys.stderr.write("Max number of cores is 20\n")
  raise SystemExit(1)

# number of runs
runs = int(sys.argv[3])

# execute runs
for i in range(0,runs):
  os.system("python run.py " + rstring + " " + pstring)
