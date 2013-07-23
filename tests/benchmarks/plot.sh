#!/bin/bash

#launch with the name of test to generate the .eps for

cat <<EOF
set terminal postscript eps color enh
set key left bottom
set output "$1.eps"
#set log y
set xlabel 'revision'
set title 'benchmark $1 - rev $2'
EOF
#echo "set terminal x11 persist"

echo "plot \\"
n=1
while read line;
do
  n=`expr $n "+" 1`
#  echo "'datatable.$1' using 1:(int(\$$n*10.0+0.5)/10.0) title '$line' w lp,\\";
  echo "'datatable.$1' using 1:$n title '$line' w lp,\\";
done < names.$1

# this forces 0.01 to be in the yrange and ends the plot list (trailing comma above)
echo "0.01 title ''"

#echo "pause -1"

