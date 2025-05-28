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

# Given a module partition name 'A' as the first argument, find all other module
# partitions 'B' it depends on either by way of direct
#   import :B;
# statements, or because a module partition 'B' so imported itself has a statement
# of the form
#   export import :C;
#
# Call this script via
#   python3 contrib/utilities/find_module_partition_dependencies.py A
# where 'A' is the name of the module partition (which must be defined in
# one of the files in the current build directory).


from glob import glob
import networkx as nx
import re
import sys

module_partition_name = sys.argv[1]
partition_to_filename = dict()

match_imports = re.compile(r"(export )?import *: *(.*);")
match_exports = re.compile(r"(export )?module dealii *: *(.*);")


# For a given source file, read through all the lines and extract the
# ones that correspond to export or import statements. For those, add
# a link to the graph.
def add_imports_for_file(module_partition_file, G):
    f = open(module_partition_file)
    lines = f.readlines()
    f.close()

    module_partition = ""
    for line in lines:
        m = match_exports.match(line)
        if m:
            module_partition = m.group(2)
            # There can only be one 'export' statement per file, so
            # stop reading after recording the mapping from partition
            # name to file name:
            partition_to_filename[module_partition] = module_partition_file
            break
    assert (
        module_partition != ""
    ), f"File <{module_partition_file} does not seem to export a module partition!"

    for line in lines:
        m = match_imports.match(line)
        if m:
            imported_partition = m.group(2)
            G.add_edge(module_partition, imported_partition)


# Create a list of all source files in the build folder
filelist = glob("module/*_module_partitions/**/*.cc*", recursive=True)
assert filelist, "Please call the script from the top-level build directory."


# For each header file, add the imports as the edges of a directed graph.
G = nx.DiGraph()
G.add_node(module_partition_name)  # Make sure the graph has at least one node
for module_partition_file in filelist:
    add_imports_for_file(module_partition_file, G)

assert (
    module_partition_name in partition_to_filename
), "The module partition given as argument could not be found in the module units found"


# Now find everything that is upstream of the module partition
# specified on the command line:
print(f"The dependencies of module partition '{module_partition_name}' are:")
for node in nx.dfs_postorder_nodes(G, source=module_partition_name):
    if node != module_partition_name:
        if node in partition_to_filename:
            print(f"  {node}, implemented in {partition_to_filename[node]}")
        else:
            print(f"  {node}, implemented in unknown file")

n_dependencies = (
    sum(1 for _ in nx.dfs_postorder_nodes(G, source=module_partition_name)) - 1
)
print(
    f"{n_dependencies} partitions are upstream dependencies for '{module_partition_name}'"
)
