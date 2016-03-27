# function to get x and y arrays from input file
def extractPlotData(input_filename):
  # get lines from input file
  input_file = open(input_filename, 'r')
  file_lines = input_file.read().split('\n')
  input_file.close()

  # initialize x and y to empty arrays
  x = []
  y = []

  # loop over lines in input file
  for file_line in file_lines:
    # ignore blank lines
    if (file_line != ''):
      # ignore lines containing "#"
      if ("#" not in file_line):
        # append to x and y arrays
        x_and_y = file_line.split()
        x.append(float(x_and_y[0]))
        y.append(float(x_and_y[1]))

  return (x,y)
