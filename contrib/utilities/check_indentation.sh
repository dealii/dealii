#!/bin/sh

if [ "${TRAVIS_PULL_REQUEST}" = "false" ]; then 
	echo "Running indentation test on master merge."
else 
	echo "Running indentation test on Pull Request #${TRAVIS_PULL_REQUEST}"
fi

export PATH=`pwd`/astyle/build/gcc/bin:$PATH

./contrib/utilities/indent
git diff
git diff-files --quiet 
