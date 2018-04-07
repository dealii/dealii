#!/bin/sh
## ---------------------------------------------------------------------
##
## Copyright (C) 2018 by the deal.II authors
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
# This script downloads and installs astyle-2.04 in programs/astyle/build/gcc.
# The installed binary is used in the 'indent' script in case astyle is
# installed by this script.
#
# This script only works on Linux, and also puts the 'astyle'
# executable into a very specific place. This is not a general purpose
# solution. A better solution is to essentially execute the commands
# below in a shell and make sure that the executable is installed into
# a central directory -- say, /usr/bin or /usr/local/bin.
#

PRG=$PWD/programs

if [ ! -d "$PRG" ]
then
    echo "create folder $PRG"
    mkdir "$PRG"
fi

# astyle
if [ ! -d "$PRG/astyle" ]
then
    echo "Downloading and installing astyle."
    mkdir "$PRG/astyle"
    wget https://downloads.sourceforge.net/project/astyle/astyle/astyle%202.04/astyle_2.04_linux.tar.gz  > /dev/null
    if echo "70b37f4853c418d1e2632612967eebf1bdb93dfbe558c51d7d013c9b4e116b60 astyle_2.04_linux.tar.gz" | sha256sum -c; then
      tar xfz astyle_2.04_linux.tar.gz -C "$PRG" > /dev/null
      cd "$PRG/astyle/build/gcc" || exit 1
      make -j4 > /dev/null
    else
      echo "*** The downloaded file has the wrong SHA256 checksum!"
      exit 1
    fi
fi
