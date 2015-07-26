#!/bin/sh

if [ "${TRAVIS_PULL_REQUEST}" = "false" ]; then 
	echo "Running indentation test on master merge."
else 
	echo "Running indentation test on Pull Request #${TRAVIS_PULL_REQUEST}"
fi

./contrib/utilities/indent
git diff
git diff-files --quiet 
