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

if [ ! -f index.html ]
then
  echo "Please run this script in the doc output directory (install/doc)"
  exit 1
fi

mkdir images >/dev/null

echo "Downloading images (press ctrl-c to cancel) ..."
cd images
{
trap "echo \"(skipping)\"" SIGINT
wget -q -nd -A png,gif -m -l 1 -np  http://www.dealii.org/images/steps/developer/
}
cd ..

echo "Patching html files ..."
sed -i 's#"http://www.dealii.org/images/steps/developer/\(step-.*\)"#"images/\1"#g' step_*.html

echo "all done!"
