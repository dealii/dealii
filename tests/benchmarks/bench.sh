#!/bin/bash

TESTS="step-22 tablehandler"

PREVREVISION="`svn info deal.II | grep Revision | sed s/Revision://`"
HEADREVISION="`svn info http://www.dealii.org/svn/dealii | grep Revision | sed s/Revision://`"
MAKECMD="make -j10"

echo "previous $PREVREVISION"
echo "HEAD: $HEADREVISION"

while [ $PREVREVISION -lt $HEADREVISION ] ; do

  NEXTREVISION=`expr $PREVREVISION "+" 100`
  echo "Updating from $PREVREVISION to $NEXTREVISION"
  cd deal.II
  svn up -r$NEXTREVISION
  if test -z "`svn diff -r$PREVREVISION:$NEXTREVISION .`" ; then
      echo "Skipping revision $NEXTREVISION" ;
      continue ;
  fi

  echo "configure"
  ./configure --disable-threads --with-petsc=no >/dev/null
  echo "compiling" 
  nice make optimized -j 10>/dev/null
  
  cd ..

  for test in $TESTS ; do
      cd $test
      echo "** working on $test"
      make clean
      make run | grep "|" > temp.txt
      cat temp.txt >>datatable.$test
      #./your_code $NEXTREVISION >>datatable.$test
      # collect info
      cd ..      
  
  done  

done

echo "DONE WITH REGRESSION TESTS ON `date`"

