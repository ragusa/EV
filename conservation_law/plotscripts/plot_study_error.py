import sys     # for command-line arguments and exit()
import os.path # for checking the existence of files
import matplotlib.pyplot as plt
from matplotlib import rc # for rendering TeX in plots
import numpy as np # for arange()

# function to extract L1 and L2 error
def extractErrorPoints(filename):
  # get lines from input file
  input_file = open(filename, 'r')
  file_lines = input_file.read().split('\n')
  input_file.close()
  # loop over lines in input file
  for file_line in file_lines:
    # ignore blank lines
    if (file_line != ''):
      # ignore lines containing "#"
      if ("#" not in file_line):
        # count number of words
        words = file_line.split()
        n_words = len(words)
        # transient convergence file
        if (n_words >= 10):
          L1 = float(words[6])
          L2 = float(words[8])
        else:
          sys.exit("Error: Need to program steady-state case\n")

  return (L1,L2)

# get command-line arguments: sys.argv[0] is script name
parameter_name = sys.argv[1]
parameter_values_string = sys.argv[2]
parameter_values = parameter_values_string.split()

# output directory
output_dir = "./output/"

# form arrays of L1 and L2 errors
L1_array = []
L2_array = []
for i in range(0,9):
  filename = output_dir + "convergence" + str(i) + ".gpl"
  if (os.path.isfile(filename)):
    (L1,L2) = extractErrorPoints(filename)
    L1_array.append(L1)
    L2_array.append(L2)

# count number of runs
n_runs = len(L1_array)

# create x positions of groups and width of bars
ind = np.arange(n_runs)
width = 0.35

# create bars
fig, ax = plt.subplots()
rects1 = ax.bar(ind,         L1_array, width, color='r')
rects2 = ax.bar(ind + width, L2_array, width, color='b')

# add some text for labels and axes ticks
ax.set_ylabel('Error')
ax.set_xticks(ind + width)
xticklabels = []
for i in range(0,n_runs):
  xticklabels.append("$" + parameter_name + " = " + parameter_values[i] + "$")
ax.set_xticklabels(xticklabels)

# legend
ax.legend((rects1[0], rects2[0]), ('$L_1$', '$L_2$'), loc='best', frameon=False)

# save figure to file
filename = "errors.pdf"
plt.savefig(filename)
