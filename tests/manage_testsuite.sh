## ---------------------------------------------------------------------
##
## Copyright (C) 2015 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

#
# Helper script used for testsuite management to configure, prune,
# regenerate and clean testsuite subprojects
#
# Usage:
#   manage_testsuite.sh configure <category> <generator> <options>
#   manage_testsuite.sh prune <category> <generator> <options>
#   manage_testsuite.sh regenerate <category> <generator> <options>
#   manage_testsuite.sh clean <category> <generator> <options>
#

set -u

CMAKE_COMMAND=cmake

CATEGORY="$2"
CATEGORY_DIR="$3"
GENERATOR="$4"
OPTIONS="$5"

configure(){
  mkdir -p "${CATEGORY}"
  cd "${CATEGORY}" && ${CMAKE_COMMAND} -G"${GENERATOR}" ${OPTIONS} "${CATEGORY_DIR}" > /dev/null 2> /dev/null
}

prune(){
  rm -rf "${PWD}/${CATEGORY}"
}

regenerate(){
  test ! -d ${CATEGORY} || ${CMAKE_COMMAND} --build "${CATEGORY}" --target regenerate > /dev/null 2> /dev/null
}

clean(){
  test ! -d ${CATEGORY} || ${CMAKE_COMMAND} --build "${CATEGORY}" --target clean > /dev/null 2> /dev/null
}

case $1 in
  configure)
    configure;;
  prune)
    prune;;
  regenerate)
    regenerate;;
  clean)
    clean;;
  *)
    exit 1;;
esac
