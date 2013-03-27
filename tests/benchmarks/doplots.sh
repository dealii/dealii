#!/bin/bash

source testlist.sh

echo "generating images..."

for test in $TESTS ; do
    LASTREV=`tail -n 1 datatable.$test | cut -f 1 -d ' '`
    ./plot.sh $test $LASTREV >script
    gnuplot script
    rm -rf script
    convert -density 150 $test.eps $test.png
done


./baselineplot.sh > script
gnuplot script
rm -rf script
convert -density 150 baseline.eps baseline.png
python interactive.py >index.html