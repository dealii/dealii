#!/bin/sh

if [ "${TRAVIS_PULL_REQUEST}" = "false" ]; then 
	echo "Running indentation test on master merge."
else 
	echo "Running indentation test on Pull Request #${TRAVIS_PULL_REQUEST}"
fi

wget http://downloads.sourceforge.net/project/astyle/astyle/astyle%202.04/astyle_2.04_linux.tar.gz > /dev/null
tar xvfz astyle_2.04_linux.tar.gz > /dev/null
cd astyle/build/gcc
make -j4 > /dev/null
cd ../../../
export PATH=`pwd`/astyle/build/gcc/bin:$PATH
./contrib/utilities/indent
git diff
git diff-files --quiet 
