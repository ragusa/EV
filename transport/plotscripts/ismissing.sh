#!/bin/bash
# helper script used by gnuplot script to print 1 if the
# file is missing and 0 if it exists
test -e ../output/$1
echo $?
