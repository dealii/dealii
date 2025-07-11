#!/usr/bin/python3

## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

# Given a C++ header file as input, read through it and output an
# equivalent interface module unit that exports a module partition
# with the declarations that are part of the header file, and that
# uses 'import' statements in place of existing '#include' directives.
#
# Call this script via
#   python3 contrib/utilities/convert_header_file_to_interface_module_unit.py in.h out.ccm


import sys
import re

import convert_to_module_units_common as header_to_partition_maps


header_file = sys.argv[1]
interface_module_unit_file = sys.argv[2]

assert header_file, "No input file name given."
assert interface_module_unit_file, "No output file name given."

m = re.search(r".*/include/deal.II/(.*).h", header_file)
interface_module_partition_name = "interface_partition_" + m.group(1).replace("/", "_")


# Read the entire header file into memory:
input = open(header_file, "r")
lines = input.readlines()
input.close()


# Regular expressions that match the preprocessor defines that we use
# to open and close the deal.II namespace. We allow for whitespace
# before and after, but nothing else in these regular
# expressions. That means that the lines in question cannot contain
# anything else -- which we will use in a few places to ensure that
# these statements are not expanded.
match_dealii_open = re.compile(r"^ *DEAL_II_NAMESPACE_OPEN *$")
match_dealii_close = re.compile(r"^ *DEAL_II_NAMESPACE_CLOSE *$")

# Then scan through it and collect the deal.II headers that are
# '#include'd, excluding "config.h" and "exception_macros.h", both of
# which declare macros and so much be #included rather than imported:
dealii_include_list = header_to_partition_maps.get_dealii_includes(lines)

# Determine which external projects that we wrap with module
# partitions are being utilized by the current file:
used_external_projects = header_to_partition_maps.get_used_external_projects(lines)


# Now open the output file and start writing into it. Start the file
# with the 'module;' text that indicates that this is going to be a
# module file. Then also include the file that defines all deal.II
# macros (including everything that's in config.h) that ensures that
# we can use preprocessor defines describing the configuration of
# deal.II. We need to explicitly include it here because we comment
# out all other #includes and so can't get at config.h transitively.
output = open(interface_module_unit_file, "w")
output.write("module;\n")
output.write("#include <deal.II/macros.h>\n\n")

# Then go through the lines of the input file again:
for line in lines:

    # If this was an '#include' directive for a deal.II header, or an
    # external project that we wrap, then just ignore/remove this
    # line.
    #
    # In source files, it is in principle acceptable to have #includes
    # because source files do not have exported partitions and so
    # processing #includes does not contribute to the bloat that
    # happens when importing interface partitions that each contain
    # dozens or hundreds of (transitive) includes. As a consequence,
    # we do allow #includes in specific circumstances, when they are
    # specifically marked with a comment. But we do not allow this
    # in interface units (header files) and so error out if someone
    # tries this.
    if header_to_partition_maps.match_dealii_includes.match(
        line
    ) or header_to_partition_maps.matches_external_header_include(line):
        if "Do not convert for module purposes" in line:
            raise "It is not allowed to escape the conversion of include statements in header files."
        else:
            pass

    # If this line contained the text DEAL_II_NAMESPACE_OPEN, then
    # prefix this text by the interface module unit start, the import
    # statements that correspond to the previously '#include'd header
    # files, and an 'export {' statement.
    #
    # Because import statements do not transitively import what the
    # imported module itself imported, we would have to annotate every
    # source file with a complete list of all module partitions (or,
    # better, corresponding header files) it needs something
    # from. This is not realistic. Rather, we emulate the semantics of
    # header files by not only importing partitions, but re-exporting
    # them as well via the 'export import :interface_partition;'
    # statement.
    elif match_dealii_open.match(line):
        output.write(
            "\nexport module dealii : " + interface_module_partition_name + ";\n\n"
        )

        for external_project in used_external_projects:
            output.write("import dealii.external." + external_project + ";\n")

        for inc in dealii_include_list:
            module_partition = "interface_partition_" + inc.replace("/", "_")
            output.write("import :" + module_partition + ";\n")

        output.write("\nexport\n{\n\n")
        output.write(line)  # Copy the previous line

    # If this line contained the text DEAL_II_NAMESPACE_CLODE, then
    # close the previous 'export {' statement:
    elif match_dealii_close.match(line):
        output.write(line)  # Copy the previous line
        output.write("\n} // close previous 'export' statement\n\n")

    # Otherwise just copy the previous line:
    else:
        output.write(line)
