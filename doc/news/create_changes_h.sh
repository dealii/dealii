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
# This script creates "changes.h.new" from the contributions in the subfolders
# of ./doc/news.
#
# The script needs to be executed as 
#   ./create_changes_h.sh
# from ./doc_news.



if test ! -d incompatibilities -o ! -d general -o ! -d specific ; then
  echo "*** This script must be run from ./doc/news!"
  exit 1
fi

OUTPUT="changes.h.new"

function process_directory
{
  echo "<ol>" >> ${OUTPUT}
  # process all entries in the right order
  ARRAY=($(ls "$1" | sort -r))
  if ! [ -z "${ARRAY[0]}" ]; then
    # treat the first item special
    # we don't want to insert a new line here
    echo -n " <li>" >> ${OUTPUT}
    # indent lines by one blank space
    sed 's/^[ ]*/ /' "$1"/"${ARRAY[0]}" >> ${OUTPUT}
    echo " </li>" >> ${OUTPUT}
    unset "ARRAY[0]"

    for f in "${ARRAY[@]}"; do
      echo "" >> ${OUTPUT}
      echo -n " <li>" >> ${OUTPUT}
      # indent lines by one blank space
      sed 's/^[ ]*/ /' "$1"/"$f" >> ${OUTPUT}
      echo " </li>" >> ${OUTPUT}
    done
  fi
  echo "</ol>" >> ${OUTPUT}
}

cat header > ${OUTPUT}

echo INCOMPATIBILITIES
cat header_incompatibilities >> ${OUTPUT}
process_directory incompatibilities

echo GENERAL
cat header_general >> ${OUTPUT}
process_directory general

echo SPECIFIC
cat header_specific >> ${OUTPUT}
process_directory specific

cat footer >> ${OUTPUT}
