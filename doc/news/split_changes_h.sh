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
# This script splits the previously used into contributions in the folders:
# incompatibilities/ general/ specific/
# The resulting files can be used in ./create_changes_h.sh to create
# changes.h anew.
#
# The script needs to be executed as 
#   ./split_changes_h.sh
# from ./doc_news.

if test ! -d incompatibilities -o ! -d specific -o ! -d general ; then
  echo "*** This script must be run from ./doc/news!"
  exit 1
fi


csplit --silent changes.h '/^<ol>\|<\/ol>$/' '{*}'

for f in `ls xx*`; do
  #remove HTML list tags
  #sed -i'' '/<ol>\|<\/ol>/!p' $f
  sed -i'' '/<ol>\|<\/ol>/d' $f
done

mv xx00 header_incompatibilities
mv xx01 incompatibilities/summary
mv xx02 header_general
mv xx03 general/summary
mv xx04 header_specific
mv xx05 specific/summary
mv xx06 footer

echo INCOMPATIBILITIES
cd incompatibilities
csplit --silent summary '/^<li>\|<\/li>$/' '{*}'
../split_summary.sh
cd ..

echo GENERAL
cd general
csplit --silent summary '/^<li>\|<\/li>$/' '{*}'
../split_summary.sh
cd ..

echo SPECIFIC
cd specific
csplit --silent summary '/^<li>\|<\/li>$/' '{*}'
../split_summary.sh
cd ..
