#!/bin/bash
## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2018 - 2024 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

#
# This script runs the clang-tidy tool on the deal.II code base.
#
#
# Usage:
# /contrib/utilities/run_clang_tidy.sh SRC_DIR OPTIONAL_CMAKE_ARGS
#   with:
#     SRC_DIR points to a deal.II source directory
#     OPTIONAL_CMAKE_ARGS are optional arguments to pass to CMake
#   make sure to run this script in an empty build directory
#
# Requirements:
# Clang 5.0.1+ and have clang, clang++, and run-clang-tidy in your path.

# grab first argument and make relative path an absolute one:
SRC=$1
SRC=$(cd "$SRC";pwd)
shift

if test ! -d "$SRC/source" -o ! -d "$SRC/include" -o ! -d "$SRC/examples" -o ! -f "$SRC/CMakeLists.txt" ; then
    echo "Usage:"
    echo "  run_clang_tidy.sh /path/to/dealII"
    exit 1
fi
echo "SRC-DIR=$SRC"

# enable MPI (to get MPI warnings)
# export compile commands (so that run-clang-tidy works)
ARGS=("-D" "DEAL_II_WITH_MPI=ON" "-D" "CMAKE_EXPORT_COMPILE_COMMANDS=ON" "-D" "CMAKE_BUILD_TYPE=Debug" "$@")

# for a list of checks, see /.clang-tidy
cat "$SRC/.clang-tidy"

if ! [ -x "$(command -v run-clang-tidy)" ] || ! [ -x "$(command -v clang++)" ]; then
    echo "make sure clang, clang++, and run-clang-tidy (part of clang) are in the path"
    exit 2
fi

CC=clang CXX=clang++ cmake "${ARGS[@]}" "$SRC" || (echo "cmake failed!"; false) || exit 2

cmake --build . --target expand_all_instantiations || (echo "make expand_all_instantiations failed!"; false) || exit 3

# generate allheaders.h
(cd include; find . -name '*.h'; cd $SRC/include/; find . -name '*.h') | grep -v allheaders.h | grep -v undefine_macros.h | sed 's|^./|#include <|' | sed 's|$|>|' >include/deal.II/allheaders.h

# finally run clang-tidy on deal.II
#
# pipe away stderr (just contains nonsensical "x warnings generated")
# pipe output to output.txt
run-clang-tidy -p . -quiet -header-filter "$SRC/include/*" -extra-arg='-DCLANG_TIDY' 2>error.txt >output.txt

# grep interesting errors and make sure we remove duplicates:
grep -E '(warning|error): ' output.txt | sort | uniq >clang-tidy.log

# if we have errors, report them and set exit status to failure
if [ -s clang-tidy.log ]; then
    cat clang-tidy.log
    exit 4
fi

echo "OK"
exit 0
