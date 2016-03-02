#!/bin/bash
## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2013 by the deal.II authors
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

# This script downloads images referenced in the tutorial steps and patches
# the URLs to point to the local files. To be run in the doc/doxygen/deal.II
# directory of an installed deal.II documentation.

if [ ! -f Tutorial.html ]
then
  echo "Please run this script in the doc output directory (<install>/doc/doxygen/deal.II)"
  exit 1
fi

mkdir -p images

echo "Downloading images (press ctrl-c to cancel) ..."
cd images
{
trap "echo \"(skipping)\"" SIGINT
wget -q -nd -A svg,png,gif -m -l 1 -np http://www.dealii.org/images/steps/developer/
}
rm robots.txt*
cd ..

echo "Patching html files ..."
sed -i 's#"http://www.dealii.org/images/steps/developer/\(step-.*\)"#"images/\1"#g' step_*.html
sed -i 's#"https://www.dealii.org/images/steps/developer/\(step-.*\)"#"images/\1"#g' step_*.html

echo "all done!"
