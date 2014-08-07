#!/bin/sh
case $1 in
#mini)
#    echo "Building and running mini tests in $2 mode."
#    mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=$2 ../ && make -j4 && make tests
#    ;;
build)
    echo "Running build tests."
    mkdir build_test 
    cd build_test
    ctest -DCMAKE_BUILD_TYPE=$2 -V -j4 -S ../tests/run_buildtest.cmake
    ;;
#tests)
#    echo "Running full testsuite."
#    mkdir full_tests
#    cd full_tests
#    ctest -DCMAKE_BUILD_TYPE=$2 -V -j4 -S ../tests/run_testsuite.cmake
#    ;;
*)
    echo "Unrecognized test type! [$1]"
    exit 1
esac

