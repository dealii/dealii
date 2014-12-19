#!/bin/sh
echo "Downloading and installing astyle."
wget http://downloads.sourceforge.net/project/astyle/astyle/astyle%202.04/astyle_2.04_linux.tar.gz  > /dev/null
tar xvfz astyle_2.04_linux.tar.gz > /dev/null
cd astyle/build/gcc
make -j4 > /dev/null
cd ../../../
export PATH=`pwd`/astyle/build/gcc/bin:$PATH
