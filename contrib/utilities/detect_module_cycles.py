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

# Read through a list of the C++ files and see if we have any cycles
# in the 'export module' or 'module' and 'import' statements of these
# files. The script here assumes that all of the files export module
# partitions of the overall module 'dealii', and import module
# partitions from the same module.
#
# Call this script via
#   python3 contrib/utilities/detect_module_cycles.py
# from the top-level directory of a build directory.

from glob import glob
import networkx as nx
import re


match_imports = re.compile(r"import *: *(.*);")
match_exports = re.compile(r"export module dealii *: *(.*);")

# For a given .cc file, read through all the lines and extract the
# ones that correspond to export or import statements. For those, add
# a link to the graph.
def add_imports_for_file(interface_module_partition_file, G):
    f = open(interface_module_partition_file)
    lines = f.readlines()
    f.close()

    for line in lines:
        m = match_exports.match(line)
        if m:
            module_partition = m.group(1)
            # There can only be one 'export' statement per file, so
            # stop reading:
            break

    for line in lines:
        m = match_imports.match(line)
        if m:
            imported_partition = m.group(1)
            G.add_edge(module_partition, imported_partition)


# Create a list of all source files in the build folder
filelist = glob("module/interface_module_partitions/**/*.cc*", recursive=True)
assert filelist, "Please call the script from the top-level directory."

# For each header file, add the imports as the edges of a directed graph.
G = nx.DiGraph()
for interface_module_partition_file in filelist:
    add_imports_for_file(interface_module_partition_file, G)

# Then figure out whether there are cycles and if so, print them:
cycles = nx.simple_cycles(G)
cycles_as_list = list(cycles)
if len(cycles_as_list) > 0:
    print(f"Cycles in the module partition graph detected!")
    for cycle in cycles_as_list:
        print(cycle)
    exit(1)
