## ---------------------------------------------------------------------
##
## Copyright (C) 2010 - 2021 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

# Remove all comments from the source file

# The output of this script is used as the source for the "The plain program"
# block at the end of each tutorial.

while (<>) {
    # Eliminate comment lines starting with "//"
    next if (m!^\s*//!);

    # Remove "//" if it is present at the end of a line. These comments exist
    # to force a specific formatting of the line for clang-format
    s!//$!!g;

    # Otherwise print the line
    print;
}	
