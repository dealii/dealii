#!/bin/bash
## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2015 - 2024 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

set -u

# Purpose: Update the copyright year of every file based on the last
#          modification recorded in the git logs
#

if test ! -d source -o ! -d include -o ! -d examples ; then
  echo "*** This script must be run from the top-level directory of deal.II."
  exit
fi

processes=1
accurate_first_year=false
until [[ "$@" == "" ]]; do
  case $1 in
    --pedantic)
      accurate_first_year=true
      shift;;
    -j)
      shift
      if [[ "$@" == "" ]]; then
        echo "Error: »-j« must be followed by a number" > /dev/stderr
        echo "Usage: update-copyright.sh [--pedantic] [-j N]" > /dev/stderr
        exit 1
      fi
      processes="${1}"
      shift;;
    *)
      echo "Error: invalid option »$1«" > /dev/stderr
      echo "Usage: update-copyright.sh [--pedantic] [-j N]" > /dev/stderr
      exit 1;;
  esac
done

#
# A shell function that updates the copyright string for a given file $1:
#

update_copyright()
{
  file="${1}"

  if ! [ -f ${file} ]; then
    echo "Skipping ${file}: not a file"
    return
  fi

  if ! head -13 ${file} | grep -q "^.. This file is part of the deal.II library.$" ; then
    echo "Skipping ${file}: no deal.II copyright header"
    return
  fi

  #
  # Get the last year this file was modified from the git log. We don't
  # want to see patches that just updated the copyright year, thus find the
  # first commit that
  #  - does not mention both the words "update" and "copyright", as well as
  #  - "Update license headers".
  #

  last_year=`git log -n 3 --date=short --format="format:%cd %s" ${file} | \
    egrep -i -v "update.*copyright|copyright.*update|Update license headers" | \
    head -n 1 | \
    perl -p -e 's/^(\d\d\d\d)-.*/\1/g;'`

  #
  # It should not happen, that the grep removes all 3 most recent commits
  # simultaneously but if it does then run the git log command again with
  # full history:
  #

  [ -z "$last_year" ] && last_year=`git log --date=short --format="format:%cd %s" ${file} | \
      egrep -i -v "update.*copyright|copyright.*update|Update license headers" | \
      head -n 1 | \
      perl -p -e 's/^(\d\d\d\d)-.*/\1/g;'`

  if [ -z "$last_year" ]; then
    echo "Skipping ${file}: internal error: could not determine last copyright year"
    return
  fi

  #
  # Get the first year this file was modified from the actual file. This is
  # fast but might be inaccurate.
  #

  first_year=`head -n 13 ${file} | egrep 'Copyright \(C\) [0-9]{4}' | \
      perl -p -e "s/.*Copyright \(C\) (\d{4}).*/\1/g;"`

  if [ -z "$first_year" ]; then
    echo "Skipping ${file}: internal error: could not determine first copyright year"
    return
  fi

  if $accurate_first_year; then
    #
    # Get the first (plausible) year this file was modified. While each file
    # (ideally) already contains a start year, experience suggests that this
    # information is typically wildly incorrect because files (and copyright
    # headers) get copied all the time. We thus grab this information from
    # git history.
    #

    #
    # First grab the oldest commit from the file history
    #

    git_first_year=`git log --reverse --date=short --format="format:%cd %s" ${file} | \
        head -n 1 | \
        perl -p -e 's/^(\d\d\d\d)-.*/\1/g;'`
    #
    # Take the minimum of what is stated in the header and the first year
    # the file was created in git:
    #
    first_year="$(( git_first_year < first_year ? git_first_year : first_year ))"

    #
    # Then, perform a more thorough search with `--diff-filter=A` to skip
    # all but the first commit in which the file was created. We try to
    # find and follow file renames with the `--follow` toggle.
    #
    # In corner cases, however, "git log --follow" can be way too
    # aggressive. As a sanity check let's restrict the admissible date
    # range to the date present in the file header, ${first_year}:
    #
    git_first_year=`git log --since="$first_year" --follow --diff-filter=A --date=short --format="format:%cd %s" ${file} | \
        tail -n 1 | \
        perl -p -e 's/^(\d\d\d\d)-.*/\1/g;'`

    #
    # If the above git command produced an output, use it. Otherwise fall
    # back to ${first_year}:
    #
    first_year="${git_first_year:-${first_year}}"
  fi

  if [ -z "$first_year" ]; then
    echo "Skipping ${file}: internal error: could not determine first copyright year"
    return
  fi

  #
  # Print a status message and update copyright line:
  #

  if [ "${first_year}" = "${last_year}" ]; then
    echo "Processing ${file}: ${last_year}"
    perl -pi -e "s/(^.. Copyright \(C\)) \d{4}( - \d{4})?(, \d{4}( - \d{4})?)*/\1 ${last_year}/g if 1..13;" ${file}
  else
    echo "Processing ${file}: ${first_year} - ${last_year}"
    perl -pi -e "s/(^.. Copyright \(C\)) \d{4}( - \d{4})?(, \d{4}( - \d{4})?)*/\1 ${first_year} - ${last_year}/g if 1..13;" ${file}
  fi
}

#
# Run copyright update in parallel:
#

process()
{
  i=0
  find ${1} -type f -regextype egrep -regex "${2}" | while read file; do
    (( i=i%processes )); (( i++==0 )) && wait
    update_copyright "${file}" &
  done
}

process "."  "CMakeLists.txt|CTestConfig.cmake" update_copyright
process "cmake contrib doc examples include source tests" ".*" update_copyright
