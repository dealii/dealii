## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2014 by the deal.II authors
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

MESSAGE(STATUS "Setting up bundled features")


IF(FEATURE_BOOST_BUNDLED_CONFIGURED)
  INSTALL(DIRECTORY ${BOOST_FOLDER}/include/boost
    DESTINATION ${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled
    COMPONENT library
    PATTERN ".svn" EXCLUDE
    )

  ADD_SUBDIRECTORY(${BOOST_FOLDER}/libs/serialization/src)
  ADD_SUBDIRECTORY(${BOOST_FOLDER}/libs/system/src)

  IF(DEAL_II_WITH_ZLIB)
    ADD_SUBDIRECTORY(${BOOST_FOLDER}/libs/iostreams/src)
  ELSE()
    MESSAGE(STATUS "BOOST::Iostreams will not be available because zlib or its header files could not be found")
  ENDIF()

  IF(DEAL_II_WITH_THREADS AND NOT DEAL_II_WITH_CXX11)
    ADD_SUBDIRECTORY(${BOOST_FOLDER}/libs/thread/src)
  ENDIF()
ENDIF()


IF(FEATURE_THREADS_BUNDLED_CONFIGURED)
  INSTALL(DIRECTORY ${TBB_FOLDER}/include/tbb
    DESTINATION ${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled
    COMPONENT library
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )

  ADD_SUBDIRECTORY(${TBB_FOLDER}/src)
ENDIF()


IF(FEATURE_UMFPACK_BUNDLED_CONFIGURED)
  ADD_SUBDIRECTORY(${UMFPACK_FOLDER}/UMFPACK/Source)
  ADD_SUBDIRECTORY(${UMFPACK_FOLDER}/AMD/Source)
ENDIF()


IF(FEATURE_MUPARSER_BUNDLED_CONFIGURED)
  ADD_SUBDIRECTORY(${MUPARSER_FOLDER})
ENDIF()



MESSAGE(STATUS "Setting up bundled features - Done")
