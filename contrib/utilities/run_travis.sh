#!/bin/sh
case $1 in
build)
    if [ "${TRAVIS_PULL_REQUEST}" = "false" ]; then 
	echo "Running build tests."
	mkdir build_test 
	cd build_test
	ctest -DCMAKE_BUILD_TYPE=$2 -V -j4 -S ../tests/run_buildtest.cmake
    else 
	echo "Build test is only run when merging to master branch. Exiting."
    fi
    ;;
indent)
     if [ "${TRAVIS_PULL_REQUEST}" = "false" ]; then 
	echo "Running indentation test on master merge."
     else 
        echo "Running indentation test on Pull Request #${TRAVIS_PULL_REQUEST}"
     fi
     wget http://downloads.sourceforge.net/project/astyle/astyle/astyle%202.04/astyle_2.04_linux.tar.gz > /dev/null
     tar xvfz astyle_2.04_linux.tar.gz > /dev/null
     pushd astyle/build/gcc
     make -j4 > /dev/null
     popd
     export PATH=`pwd`/astyle/build/gcc/bin:$PATH
     ./contrib/utilities/indent
     git diff-files --quiet || (git diff && failing_missing_command)	
     ;;
*)
    echo "Unrecognized test type! [$1]"
    exit 1
esac

