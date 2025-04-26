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
            count = count - 1
        if count < 0:
            sys.exit("Error in file '%s' in line %d" % (filename, lineno))
        lineno = lineno + 1

    if count != 0:
        sys.exit("Error: missing closing braces in file '%s'" % (filename))

    return


# have empty lines in html tables?
def check_empty_lines_in_tables(lines):
    count = 0
    lineno = 1
    for l in lines:
        if "<table" in l:
            count = count + 1
        elif "</table>" in l:
            count = count - 1
        if count < 0:
            sys.exit("Error in file '%s' in line %d" % (filename, lineno))

        if count == 1:
            if l.strip() == "":
                sys.exit(
                    "Error: empty line inside html table in file '%s' in line %d"
                    % (filename, lineno)
                )

        lineno = lineno + 1

    if count != 0:
        sys.exit("Error: mismatched html table tags in file '%s'" % (filename))

    return


# have more than one header with the same name in a tutorial?
def check_multiple_defined_headers(lines):
    headers = []
    for l in lines:
        # this may break if a header splits across a line
        if "<h" in l:
            # convert thisheader to the doxygen header name equivalent
            thisheader = l  # take this line
            thisheader = thisheader.replace(" ", "")  # remove all whitespace
            thisheader = thisheader.split("<h")[1]  # remove header tag '<h'

            if thisheader[0] == "2":
                # We do not use <h2> headers in tutorials because our script
                # does not put them in the table of contents (for historical
                # reasons). Instead, please use <h3>, <h4>, etc..
                sys.exit(
                    "Error: Header <h2> detected in file '%s'. This is"
                    " not allowed in tutorial programs. Please use <h3>"
                    " instead." % (filename)
                )

            if thisheader[0] in [
                "1",
                "2",
                "3",
                "4",
                "5",
            ]:  # make sure the next character is 1-5
                thisheader = thisheader[2:]  # remove '*>'
                if (
                    len(thisheader.split("</h")) == 2
                ):  # check for close header on the same line
                    thisheader = thisheader.split("</h")[0]  # remove close header tag
                    thisheader = "".join(
                        ch for ch in thisheader if ch.isalnum()
                    )  # remove nonalphanumeric
                    if thisheader in headers:
                        sys.exit(
                            "Error: repeated header title '%s' in file group '%s'"
                            % (thisheader, filename)
                        )
                    else:
                        headers.append(thisheader)
                else:
                    sys.exit(
                        "Error: mismatched html header flags in file group '%s'"
                        % (filename)
                    )
            # do not check else here because it may match with template arguments like '<hsize_t>'

        # detect @sect[1-5]{ ... }
        if "@sect" in l:
            thisheader = l  # take this line
            thisheader = thisheader.replace(" ", "")  # remove all whitespace
            thisheader = thisheader.split("@sect")[1]  # remove '@sect'
            if thisheader[0] in [
                "1",
                "2",
                "3",
                "4",
                "5",
            ]:  # make sure the next character is 1-5
                thisheader = thisheader[2:]  # remove the number
                thisheader = "".join(
                    ch for ch in thisheader if ch.isalnum()
                )  # remove nonalphanumeric

                if thisheader in headers:
                    sys.exit(
                        "Error: repeated header title '%s' in file group '%s'"
                        % (thisheader, filename)
                    )
                else:
                    headers.append(thisheader)

    return


args = sys.argv
if len(args) != 2:
    print("Usage: %s <filename>" % (sys.argv[0]))
    sys.exit("Error: please pass exactly one filename to this script")

filename = args[1]

f = open(filename)
lines = f.readlines()
f.close()

check_doxygen_groups(lines)
check_empty_lines_in_tables(lines)

# if it's intro.dox, combine it with step-**.cc and results.dox and check headers
if "intro.dox" in filename:
    # get stepname.cc
    stepname = filename.split("/doc", 1)[0]
    stepname = stepname + "/" + stepname.split("/", 1)[1] + ".cc"

    # get results.dox
    resultsname = filename.split("intro.dox", 1)[0] + "results.dox"

    # open the '.cc' file:
    fstep = open(stepname)

    lines += fstep.readlines()
    fstep.close()

    fres = open(resultsname)
    lines += fres.readlines()
    fres.close()

    # check the file group
    check_multiple_defined_headers(lines)
