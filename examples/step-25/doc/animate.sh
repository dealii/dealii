#!/bin/sh
## ---------------------------------------------------------------------
##
## Copyright (C) 2006 - 2013 by the deal.II authors
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
# Author: Ivan Christov, Wolfgang Bangerth, Texas A&M University, 2006
#

dim=1
tstep=0
tstepinc=20
numtsteps=6300
batchfile='animation.plt'

# optional clean up...all these files will get overwritten
rm "solution-"$dim"d-"*".png"
rm "solution-"$dim"d.gif"

# always gotta delete this one, though
rm $batchfile

# generate the gnuplot batch script to plot all the desire time steps
while [ $tstep -lt $numtsteps ]
do
    if [ $tstep -lt 10 ]; then
	ztstep="0000"$tstep
    elif [ $tstep -lt 100 ]; then
	ztstep="000"$tstep
    elif [ $tstep -lt 1000 ]; then
	ztstep="00"$tstep
    elif [ $tstep -lt 10000 ]; then
	ztstep="0"$tstep
    else
	ztstep=$tstep
    fi

    echo "call \"plot.plt\" $dim $ztstep" >> $batchfile
    
    let tstep=tstep+$tstepinc
done

gnuplot -persist $batchfile

# use ImageMagick to create an animated gif from the PNG files
convert -delay 0 -loop 0 "solution-"$dim"d-*.png" "solution-"$dim"d.gif"
