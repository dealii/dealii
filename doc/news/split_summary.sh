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
# This script is used in ./split_changes_h.sh to split entries 
# in the respective into small files.
#
# The script needs to be executed as 
#   ../split_summary.sh
# from one of the subfolders of doc/news that contains a 'summary' file
# created by ./split_changes.h.
#

if test ! -f summary ; then
  echo "*** No 'summary' file found!"
  exit 1
fi


rm summary
for f in `ls xx*`; do
  #remove HTML list tags
  sed -i'' 's/<li>\|<\/li>\|<ol>\|<\/ol>//g' $f
  #remove trailing whitespace
  sed -i'' 's/^[ \t]*//' $f
  #remove empty lines
  sed -i'' '/^\s*$/d' $f
  #only consider non-empty files 
  if [[ -s $f ]] ; then
    cat $f > tmp
    DATE=`sed -n -r 's/^.*([0-9]{4})\/([0-1][0-9])\/([0-3][0-9]).*$/\1\2\3/p' $f`
    TMP=`tail $f | sed -n -r 's/.*\(([A-Za-z0-9, ]*).*[0-9]{4}\/[0-1][0-9]\/([0-3][0-9])\).*/\1/p'`
    NAME=`echo $TMP | sed -r 's/[ ,]+//g'`
    OLDNAME=$NAME
    COUNTER=0
    while [[ -s ${DATE}${NAME} ]] ; do
      COUNTER=$((COUNTER+1)) 
      NAME="${OLDNAME}_${COUNTER}"
    done
    mv tmp ${DATE}${NAME}
  fi
  rm $f
done
