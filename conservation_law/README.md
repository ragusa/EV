# Viewing this README in web browser

This README can be viewed in a web browser using the script `./bin/md`:
```
./bin/md README.md
```
This script assumes that the following are available terminal commands:

* `pandoc`, for converting markdown syntax to html, and
* `google-chrome`, for viewing the resulting html.

# Directory structure

The directory structure is as follows:

Directory       | Contents
----------------|-----------------------------
`./bin`         | Useful scripts
`./CMakeFiles`  | `CMake` output files
`./doxygen`     | `doxygen` output files
`./include`     | Header files
`./input`       | Run input files
`./mesh`        | Mesh input files
`./output`      | Output files
`./plotscripts` | Plot scripts
`./problems`    | Problem input files
`./src`         | Source files
`./studies`     | Study input files
`./templates`   | Template files (for creating a new class)
`./Testing`     | `Ctest` output files
`./tests`       | unit/regression test input and output files

There are currently 2 main files, one for explicit runs: `./main.cc`, and the
other for implicit transport `./implicit_transport.cc`.

# Running executables

1. Run `cmake`:
```
cmake .
```
2. Run `make`:
```
make
```
3. Edit the input file corresponding to the conservation law for which you
   wish to run, e.g., `./input/input_burgers.prm` or `./input/input_euler.prm`.
4. Run the executable created by make, e.g., `burgers` or `euler`:
```
./burgers
```

# Plotting results

Run the script `./bin/plot`:
```
./bin/plot
```

This script calls the script `./bin/plot_explicit` with arguments inferred by the
contents of `CMakeLists.txt` and the appropriate input file. To call
`plot_explicit` directly, call it without arguments to see the correct usage.

# Running unit/regression tests

Currently all tests are 1-D, so one must first check that the dimension constant
in `main.cc` and `implicit_transport.cc` is equal to 1. Then, tests can be run
with the following command:
```
ctest
```
Note there is no need to make sure that executables are compiled; `ctest` will
do this automatically.

# Viewing doxygen

## Using the script

The script `./bin/dox` executes the above and switches focus to the web browser.
It assumes that Google Chrome (`google-chrome`) and `wmctrl` are installed:
```
./bin/dox
```

## Manually viewing doxygen

Doxygen documentation can be compiled as follows:
```
doxygen Doxyfile
```
Then `./doxygen/html/index.html` can be opened in your favorite web browser.

