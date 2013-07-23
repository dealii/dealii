#!/bin/bash

source testlist.sh

PREVREVISION="`svn info deal.II | grep Revision | sed s/Revision://`"
MAKECMD="nice make -j10 install"
#MAKECMD="nice make -j10 optimized"
export MAKECMD

echo "testing $PREVREVISION"

  echo "configure"
  cd build
  cmake ../deal.II || exit 2
  #./configure --disable-threads --with-petsc=no || exit 2
  echo "compiling" 
  $MAKECMD || exit 3
  
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
      rm -rf names.test
      fi ;
      ./../gettimes/gettimes $PREVREVISION >>../datatable.$test
      cd ..      
  
  done  
