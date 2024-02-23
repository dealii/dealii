## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2010 - 2023 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

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
