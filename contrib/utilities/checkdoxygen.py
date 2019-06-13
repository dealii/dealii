#!/usr/bin/python

# run this script on all headers of deal.II to check various doxygen related problems

# for example:
# find doc examples include \( -name "*.h" -o -name "*.dox" \) -print | xargs -n 1 contrib/utilities/checkdoxygen.py


import sys

# have matching @{ @} pairs? (doxygen groups)
def check_doxygen_groups(lines):
    count = 0
    lineno = 1
    for l in lines:
        if "@{" in l:
            count = count + 1
        elif "@}" in l:
            count = count -1
        if count < 0:
            sys.exit("Error in file '%s' in line %d"%(args[0], lineno))
        lineno = lineno + 1

    if count != 0:
        sys.exit("Error: missing closing braces in file '%s'"%(args[0]))

    return

# have empty lines in html tables?
def check_empty_lines_in_tables(lines):
    count = 0
    lineno = 1
    for l in lines:
        if "<table" in l:
            count = count + 1
        elif "</table>" in l:
            count = count -1
        if count < 0:
            sys.exit("Error in file '%s' in line %d"%(args[0], lineno))

        if count == 1:
            if l.strip() == "":
                sys.exit("Error: empty line inside html table in file '%s' in line %d"%(args[0], lineno))

        lineno = lineno + 1

    if count != 0:
        sys.exit("Error: mismatched html table tags in file '%s'"%(args[0]))

    return


args = sys.argv
args.pop(0)

f = open(args[0])
lines = f.readlines()
f.close()

check_doxygen_groups(lines)
check_empty_lines_in_tables(lines)


