#!/bin/bash

TESTS=step-22

PREVREVISION="`svn info deal.II | grep Revision | sed s/Revision://`"
HEADREVISION="`svn info http://www.dealii.org/svn/dealii | grep Revision | sed s/Revision://`"
MAKECMD="make -j16"

echo "previous $PREVREVISION"
echo "HEAD: $HEADREVISION"

while [ $PREVREVISION -lt $HEADREVISION ] ; do

  NEXTREVISION=`expr $PREVREVISION "+" 1`
  echo "Updating from $PREVREVISION to $NEXTREVISION"
  pause
  cd deal.II
  svn up deal.II -r$NEXTREVISION
  echo "configure"
  ./configure --disable-threads --with-petsc=no >/dev/null
  echo "compiling" 
  nice make optimized -j 10>/dev/null
  
  $MAKECMD optimized
  cd ..

  for test in $TESTS ; do
      cd $test
      echo "** working on $test"

      make run
      # collect info
      cd ..      
  
  done  

fi

echo "DONE WITH REGRESSION TESTS ON `date`"

