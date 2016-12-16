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
# This script splits the previously used "changes.h" into contributions
# in the folders "incompatibilities", "major" and "minor".
# The resulting files can be used in ./create_changes_h.sh to create
# changes.h anew.
#
# The script needs to be executed as 
#   ./split_changes_h.sh
# from ./doc/news/changes.

if test ! -d incompatibilities -o ! -d minor -o ! -d major ; then
  echo "*** This script must be run from ./doc/news/changes!"
  exit 1
fi

if test ! -f ../changes.h ; then
  echo "*** The file '../changes.h' does not exist!"
  exit 1
fi



csplit --silent ../changes.h '/^<ol>\|<\/ol>$/' '{*}'

for f in xx*; do
  #remove HTML list tags
  sed -i'' '/<ol>\|<\/ol>/d' "$f"
done

mv xx00 header_incompatibilities
mv xx01 incompatibilities/summary
mv xx02 header_major
mv xx03 major/summary
mv xx04 header_minor
mv xx05 minor/summary
mv xx06 footer

csplit --silent header_incompatibilities '/^<!--.*$/' '{*}'
mv xx00 header
mv xx01 header_incompatibilities

echo INCOMPATIBILITIES
cd incompatibilities || exit
csplit --silent summary '/^<li>\|<\/li>$/' '{*}'
../split_summary.sh
cd ..

echo GENERAL
cd major || exit
csplit --silent summary '/^<li>\|<\/li>$/' '{*}'
../split_summary.sh
cd ..

echo SPECIFIC
cd minor || exit
csplit --silent summary '/^<li>\|<\/li>$/' '{*}'
../split_summary.sh
cd ..
