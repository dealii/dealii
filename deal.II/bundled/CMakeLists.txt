## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2012 - 2013 by the deal.II authors
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
# Compile and install enabled bundled features:
#

MESSAGE(STATUS "Setup bundled features")


IF(FEATURE_BOOST_BUNDLED_CONFIGURED)
  INSTALL(DIRECTORY ${BOOST_FOLDER}/include/boost
    DESTINATION ${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled
    COMPONENT library
    PATTERN ".svn" EXCLUDE
    )

  FILE(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/deal_ii_source_includes
    "${BOOST_FOLDER}/include\n"
    )

  ADD_SUBDIRECTORY(${BOOST_FOLDER}/libs/serialization/src)

  IF(DEAL_II_WITH_THREADS AND NOT DEAL_II_USE_CXX11)
    #
    # If the C++ compiler doesn't completely support the C++11 standard
    # (and consequently we can't use std::thread, std::mutex, etc), then
    # include all the files that form BOOST's thread implementation so that
    # we don't have to build BOOST itself only to get at this small part of
    # it. it also ensures that we use the correct compiler and flags
    #
    ADD_SUBDIRECTORY(${BOOST_FOLDER}/libs/thread/src)
  ENDIF()
ENDIF()


IF(FEATURE_FUNCTIONPARSER_BUNDLED_CONFIGURED)
  ADD_SUBDIRECTORY(${FUNCTIONPARSER_FOLDER})
ENDIF()


IF(FEATURE_THREADS_BUNDLED_CONFIGURED)
  INSTALL(DIRECTORY ${TBB_FOLDER}/include/tbb
    DESTINATION ${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled
    COMPONENT library
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )

  FILE(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/deal_ii_source_includes
    "${TBB_FOLDER}/include\n"
    )

  ADD_SUBDIRECTORY(${TBB_FOLDER}/src)
ENDIF()


IF(FEATURE_UMFPACK_BUNDLED_CONFIGURED)
  ADD_SUBDIRECTORY(${UMFPACK_FOLDER}/UMFPACK/Source)
  ADD_SUBDIRECTORY(${UMFPACK_FOLDER}/AMD/Source)
ENDIF()


MESSAGE(STATUS "Setup bundled features - Done")
