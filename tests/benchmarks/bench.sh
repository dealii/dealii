#!/bin/bash

TESTS="step-22 tablehandler"

PREVREVISION="`svn info deal.II | grep Revision | sed s/Revision://`"
HEADREVISION="`svn info http://www.dealii.org/svn/dealii | grep Revision | sed s/Revision://`"
MAKECMD="make -j2"
export MAKECMD

echo "previous $PREVREVISION"
echo "HEAD: $HEADREVISION"

while [ $PREVREVISION -lt $HEADREVISION ] ; do

  NEXTREVISION=`expr $PREVREVISION "+" 100`
  echo "Updating from $PREVREVISION to $NEXTREVISION"
  cd deal.II
  svn up -r$NEXTREVISION >/dev/null
  if test -z "`svn diff -r$PREVREVISION:$NEXTREVISION .`" ; then
      echo "Skipping revision $NEXTREVISION" ;
      continue ;
  fi

  PREVREVISION=$NEXTREVISION

  echo "configure"
  ./configure --disable-threads --with-petsc=no >/dev/null
  echo "compiling" 
  $MAKECMD optimized>/dev/null
  
  cd ..

  for test in $TESTS ; do
      cd $test
      echo "** working on $test"
      make clean >/dev/null
      make run | grep "|" > temp.txt
      ./../gettimes/gettimes > names.test
      if [[ -s names.test ]] ; then
      cp names.test ../names.$test
      fi ;
      ./../gettimes/gettimes $NEXTREVISION >>../datatable.$test
      cd ..      
  
  done  

done

echo "DONE WITH REGRESSION TESTS ON `date`"
