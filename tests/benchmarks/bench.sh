#!/bin/bash

source testlist.sh

PREVREVISION="`svn info deal.II | grep Revision | sed s/Revision://`"
HEADREVISION="`svn info https://svn.dealii.org/trunk/deal.II | grep Revision | sed s/Revision://`"
MAKECMD="nice make -j20"
export MAKECMD

echo "previous $PREVREVISION"
echo "HEAD: $HEADREVISION"

while [ $PREVREVISION -lt $HEADREVISION ] ; do

  NEXTREVISION=`expr $PREVREVISION "+" 10`
  echo "Updating from $PREVREVISION to $NEXTREVISION"
  cd deal.II
  svn up -r$NEXTREVISION || exit 1
  if test -z "`svn diff -r$PREVREVISION:$NEXTREVISION .`" ; then
      echo "Skipping revision $NEXTREVISION" ;
      PREVREVISION=$NEXTREVISION
      cd ..
      continue ;
  fi

  PREVREVISION=$NEXTREVISION

  echo "configure"
  ./configure --disable-threads --with-petsc=no || exit 2
  echo "compiling" 
  $MAKECMD optimized || exit 3
  
  cd ..

  for test in $TESTS ; do
      cd $test
      echo "** working on $test"
      make clean >/dev/null
      echo -n "" > temp.txt
      for a in {1..5}; do
	  echo "*" >> temp.txt
          make run | grep "|" >> temp.txt
      done
      ./../gettimes/gettimes > names.test
      if [[ -s names.test ]] ; then
      words=`wc -w names.test | cut -f1 -d' '`
      if [ "$words" -gt "0" ] ; then
	  cp names.test ../names.$test
      fi ;
      fi ;
      ./../gettimes/gettimes $NEXTREVISION >>../datatable.$test
      cd ..      
  
  done  

done

echo "DONE WITH REGRESSION TESTS ON `date`"
