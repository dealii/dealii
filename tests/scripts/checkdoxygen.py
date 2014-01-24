#!/usr/bin/python

# run this script on all headers of deal.II to check that
# all doxygen groups have matching @{ @} pairs.
# for example:
# find include -name "*h" -print | xargs -n 1 ../tests/scripts/checkdoxygen.py

import sys

args=sys.argv
args.pop(0)

f = open(args[0])
lines = f.readlines()
f.close()


count = 0
lineno = 1
for l in lines:
    if "@{" in l:
            count = count + 1
    elif "@}" in l:
             count = count -1
             if (count < 0):
                 sys.exit("Error in file '%s' in line %d"%(args[0],lineno));
    lineno = lineno + 1

if (count != 0):
    sys.exit("Error: missing closing braces in file '%s'"%(args[0]));
