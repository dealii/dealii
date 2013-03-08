#!/bin/bash


source testlist.sh


cat <<EOF
set terminal postscript eps color enh
set key left top
set output "baseline.eps"
#set log y
set xlabel 'revision'
set ylabel 'percentage slowdown to baseline'
set title 'benchmark baselines'
EOF
#echo "set terminal x11 persist"
echo "set xrange [28000:*]"
echo "set yrange [*:300]"
echo "plot \\"

for test in $TESTS ; do
    
    col=1
    while read line;
    do
	col=`expr $col "+" 1`
	baseline=`head -n 1 datatable.$test | cut -f $col -d ' '`
	echo "'datatable.$test' using 1:((\$$col-$baseline)/$baseline*100.0) title '$test - $line' w lp,\\";
    done < names.$test

done

echo "0.0 w l title 'baseline'"

#echo "pause -1"

