#!/bin/sh
## ---------------------------------------------------------------------
##
## Copyright (C) 2014 - 2016 by the deal.II authors
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
# satisfies our indentation standards.
#
# It does so by running the 'indent' script (located in the current
# directory), calling 'git diff' to show what differences exist between
# the correctly indented code and what is in the git index (which is
# typically what is in the last commit), and then running a command
# that either returns success or failure, depending on whether or not
# there are differences. The continuous integration services return
# a failure code for a pull request if this script returns a failure.
#

if [ "${TRAVIS_PULL_REQUEST}" = "false" ]; then 
	echo "Running indentation test on master merge."
else 
	echo "Running indentation test on Pull Request #${TRAVIS_PULL_REQUEST}"
fi

./contrib/utilities/indent
git diff
git diff-files --quiet 
