#!/bin/bash
## ---------------------------------------------------------------------
##
## Copyright (C) 2016 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

#
# This script creates changes.h from the contributions in the subfolders
# of ./doc/news.
#
# The script needs to be executed as 
#   ./create_changes_h.sh
# from ./doc_news.



if test ! -d incompatibilities -o ! -d general -o ! -d specific ; then
  echo "*** This script must be run from ./doc/news!"
  exit 1
fi

echo INCOMPATIBILITIES
cat header_incompatibilities > changes.h
echo "<ol>" >> changes.h	
ARRAY=($(ls incompatibilities | sort -r))
if ! [ -z "$ARRAY" ]; then
  echo -n " <li>" >> changes.h
  sed 's/^/ /' incompatibilities/${ARRAY[0]} >> changes.h
  echo " </li>" >> changes.h
  unset ARRAY[0]

  for f in ${ARRAY[@]}; do
    echo "" >> changes.h
    echo -n " <li>" >> changes.h
    sed 's/^/ /' incompatibilities/$f >> changes.h
    echo " </li>" >> changes.h
  done
fi
echo "</ol>" >> changes.h

echo GENERAL
cat header_general >> changes.h
echo "<ol>" >> changes.h
ARRAY=($(ls general | sort -r))
if ! [ -z $ARRAY ]; then
  echo -n " <li>" >> changes.h
  sed 's/^/ /' general/${ARRAY[0]} >> changes.h
  echo " </li>" >> changes.h
  unset ARRAY[0]

  for f in ${ARRAY[@]}; do
    echo "" >> changes.h
    echo -n " <li>" >> changes.h
    sed 's/^/ /' general/$f >> changes.h
    echo " </li>" >> changes.h
  done
fi
echo "</ol>" >> changes.h

echo SPECIFIC
cat header_specific >> changes.h
echo "<ol>" >> changes.h
ARRAY=($(ls specific | sort -r))
if ! [ -z $ARRAY ]; then
  echo -n " <li>" >> changes.h
  sed 's/^/ /' specific/${ARRAY[0]} >> changes.h
  echo " </li>" >> changes.h
  unset ARRAY[0]

  for f in ${ARRAY[@]}; do
    echo "" >> changes.h
    echo -n " <li>" >> changes.h
    sed 's/^/ /' specific/$f >> changes.h
    echo " </li>" >> changes.h
  done
fi
echo "</ol>" >> changes.h

cat footer >> changes.h
