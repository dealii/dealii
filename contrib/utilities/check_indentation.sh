#!/bin/sh
## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2014 - 2022 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

#
# This is a script that is used by the continuous integration servers
# to make sure that the currently checked out version of a git repository
# satisfies our "indentation" standards. This does no longer only cover
# indentation, but other automated checks as well.
#
# WARNING: The continuous integration services return a failure code for a
# pull request if this script returns a failure, so the return value of this
# script is important.
#

# Run indent-all and fail if script fails:
./contrib/utilities/indent-all || exit $?

# Run lowercase_cmake and fail if script fails:
./contrib/utilities/lowercase_cmake || exit $?

# Show the diff in the output:
git diff

# Make this script fail if any changes were applied by indent above:
git diff-files --quiet || exit $?

# Run various checks involving doxygen documentation:
./contrib/utilities/check_doxygen.sh || exit $?

# Check for utf8 encoding:
./contrib/utilities/check_encoding.py || exit $?

# Success!
exit 0
