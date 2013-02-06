#!/bin/bash


while :
do
   ./bench.sh

   source testlist.sh

   for test in $TESTS ; do
       LASTREV=`tail -n 1 datatable.$test | cut -f 1 -d ' '`
       ./plot.sh $test $LASTREV >script
       gnuplot script
       rm -rf script
       convert -density 150 $test.eps $test.png
   done
   
   ./sync.sh

   sleep 3600
done