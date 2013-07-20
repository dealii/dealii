## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2013 by the deal.II authors
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

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

INCLUDE_DIRECTORIES(
  BEFORE
  ${CMAKE_SOURCE_DIR}/include/
  ${CMAKE_BINARY_DIR}/include/
  )

ADD_EXECUTABLE(test_affinity
  test_affinity.cc
  )

IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  SET(_build "DEBUG")
ELSE()
  SET(_build "RELEASE")
ENDIF()

SET_TARGET_PROPERTIES(test_affinity PROPERTIES
  LINK_FLAGS "${DEAL_II_LINKER_FLAGS} ${DEAL_II_LINKER_FLAGS_${_build}}"
  COMPILE_DEFINITIONS "${DEAL_II_DEFINITIONS};${DEAL_II_DEFINITIONS_${_build}}"
  COMPILE_FLAGS "${DEAL_II_CXX_FLAGS_${_build}}"
  LINKER_LANGUAGE "CXX"
  )

TARGET_LINK_LIBRARIES(test_affinity ${DEAL_II_BASE_NAME}${DEAL_II_${_build}_SUFFIX})

IF(NOT CMAKE_CROSSCOMPILING)
  ADD_CUSTOM_TARGET(run_test_affinity COMMAND test_affinity)
ENDIF()
