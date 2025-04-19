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

# Given the name of a source or header file, find all deal.II header
# files it depends on directly or indirectly.
#
# Call this script via
#   python3 contrib/utilities/find_header_dependencies.py <filename>
# where '<filename>' is the name of the source file to be analyzed.


from glob import glob
import networkx as nx
import re
import sys

source_file = sys.argv[1]

match_dealii_includes = re.compile(r"# *include *<(deal.II/.*)>")


# For a given header or source file, read through all the lines and
# extract the ones that correspond to further #include statements. For
# those, add a link to the graph.
def add_includes_for_file(file, G):
    f = open(file)
    lines = f.readlines()
    f.close()

    for line in lines:
        m = match_dealii_includes.match(line)
        if m:
            included_file = m.group(1)
            G.add_edge(file, "include/" + included_file)


# Create a list of all source files in the build folder
filelist = glob("include/deal.II/**/*.h", recursive=True)
assert filelist, "Please call the script from the top-level deal.II directory."
if source_file not in filelist:
    filelist.append(source_file)


# For each header file plus the one provided on the command line, add
# the includes as the edges of a directed graph.
G = nx.DiGraph()
G.add_node(source_file)  # Make sure the graph has at least one node
for file in filelist:
    add_includes_for_file(file, G)

# Now find everything that is upstream of the module partition
# specified on the command line:
print(f"The dependencies of module partition '{source_file}' are:")
for node in nx.dfs_postorder_nodes(G, source=source_file):
    if node != source_file:
        print(f"  {node}")

n_dependencies = sum(1 for _ in nx.dfs_postorder_nodes(G, source=source_file)) - 1
print(f"{n_dependencies} files are upstream dependencies for '{source_file}'")
