#!/bin/bash

PREVREVISION="`svn info deal.II | grep Revision | sed s/Revision://`"
HEADREVISION="`svn info https://svn.dealii.org/trunk/deal.II | grep Revision | sed s/Revision://`"

echo "previous $PREVREVISION"
echo "HEAD: $HEADREVISION"

while [ $PREVREVISION -lt $HEADREVISION ] ; do

  NEXTREVISION=`expr $PREVREVISION "+" 25`
  echo "Updating from $PREVREVISION to $NEXTREVISION"
  cd deal.II
  svn up -r$NEXTREVISION || exit 1
  if test -z "`svn diff -r$PREVREVISION:$NEXTREVISION .`" ; then
      echo "Skipping revision $NEXTREVISION" ;
      PREVREVISION=$NEXTREVISION
      cd ..
      continue ;
  fi

  cd ..      
  PREVREVISION=$NEXTREVISION
  . benchrev.sh 

done

echo "DONE WITH REGRESSION TESTS ON `date`"
