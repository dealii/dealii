#!/bin/sh
## ---------------------------------------------------------------------
##
## Copyright (C) 2014 - 2018 by the deal.II authors
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

if [ "${TRAVIS_PULL_REQUEST}" = "false" ]; then 
	echo "Running indentation test on master merge."
else 
	echo "Running indentation test on Pull Request #${TRAVIS_PULL_REQUEST}"
fi

# Run indent-all and fail if script fails:
./contrib/utilities/indent-all || exit $?

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
