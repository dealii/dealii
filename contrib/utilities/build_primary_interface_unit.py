#!/usr/bin/python3

## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2014 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

# Given a list of interface module units that each contain an exported
# module partition, create a primary module interface unit that
# creates the 'dealii' module.
#
# Call this script via
#   python3 contrib/utilities/convert_header_file_to_interface_module_unit.py <list of interface module unit files>
# The program outputs the primary module interface unit to the console.


import sys
import re


match_export = re.compile(r"^export module dealii *: *(.*);")


# Print the header of the primary module partition:
print("module;")
print("export module dealii;")

# Go through all input files and check their exported module partitions:
for module_input_file in sys.argv[1:]:
    input = open(module_input_file, "r")

    # Read through the lines of the file and see where it exports a
    # module partition. Then 'export import' that partition:
    for line in input:
        m = match_export.match(line)
        if m:
            print("export import :" + m.group(1) + ";")

            # A file can only contain a single interface module
            # partition. So once we found one, we can stop parsing.
            break
