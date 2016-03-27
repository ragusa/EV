import sys # for command-line arguments
import os.path # for checking the existence of files
import matplotlib.pyplot as plt
from matplotlib import rc # for rendering TeX in plots
from plot_utilities import extractPlotData

# get command-line arguments: sys.argv[0] is script name
parameter_name = sys.argv[1]
parameter_values_string = sys.argv[2]
parameter_values = parameter_values_string.split()

# output directory
output_dir = "./output/"

# set up plot figure
plt.figure() # create new figure
plt.rc('text', usetex=True)          # use TeX to generate text
plt.rc('font', family='sans-serif')  # use sans-serif font family
plt.xlabel('$x$')
plt.ylabel('$u(x)$')

# plot exact solution if it exists
exact_solution_filename = output_dir + "solution_exact.gpl"
if (os.path.isfile(exact_solution_filename)):
  (x_exact,y_exact) = extractPlotData(exact_solution_filename)
  plt.plot(x_exact, y_exact, 'k-', label='Exact')

# loop over possible run numbers and plot
for i in range(0,9):
  filename = output_dir + "solution" + str(i) + ".gpl"
  if (os.path.isfile(filename)):
    (x,y) = extractPlotData(filename)
    plotlabel = "$" + parameter_name + " = " + parameter_values[i] + "$"
    plt.plot(x, y, label=plotlabel)
  else:
    break

# legend
plt.legend(loc='best', frameon=False)

# save figure to file
filename = "solutions.pdf"
plt.savefig(filename)
