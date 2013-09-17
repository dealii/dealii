#!/bin/bash

# set the version we want to start with. the first revision that
# can be used with this script is r24500
REV=30125

rm -rf deal.II
svn co -r $REV https://svn.dealii.org/trunk/deal.II

mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DDEAL_II_WITH_THREADS=OFF -DCMAKE_INSTALL_PREFIX=`pwd`/../installed ../deal.II
make install -j 10
cd ..

cd gettimes
make
cd ..

source testlist.sh
for test in $TESTS ; do
      cd $test
      echo "** cmake for $test:"
      cmake -DDEAL_II_DIR=`pwd`/../installed
      make release
      cd ..
done

