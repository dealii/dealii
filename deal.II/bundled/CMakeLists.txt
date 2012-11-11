#####
##
## Copyright (C) 2012 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

#
# Export information about bundled library locations and do the actual
# setup of compilation targets and installation here:
#


#
# Boost C++ libraries
#
SET(FEATURE_BOOST_HAVE_BUNDLED TRUE)
OPTION(DEAL_II_FORCE_BUNDLED_BOOST
  "Always use the bundled boost library instead of an external one."
  OFF)
SET(BOOST_FOLDER "${CMAKE_SOURCE_DIR}/bundled/boost-1.49.0")

IF(FEATURE_BOOST_BUNDLED_CONFIGURED)
  INSTALL(DIRECTORY ${BOOST_FOLDER}/include/boost
    DESTINATION ${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled
    COMPONENT library
    PATTERN ".svn" EXCLUDE
    )

  ADD_SUBDIRECTORY(${BOOST_FOLDER}/libs/serialization/src)

  IF( DEAL_II_USE_MT AND NOT DEAL_II_CAN_USE_CXX11)
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


#
# FunctionParser project:
#
SET(FEATURE_FUNCTIONPARSER_HAVE_BUNDLED TRUE)
OPTION(DEAL_II_FORCE_BUNDLED_FUNCTIONPARSER
  "Always use the bundled functionparser library instead of an external one."
  OFF)
SET(FUNCTIONPARSER_FOLDER "${CMAKE_SOURCE_DIR}/bundled/functionparser/")

IF(FEATURE_FUNCTIONPARSER_BUNDLED_CONFIGURED)
  ADD_SUBDIRECTORY(${FUNCTIONPARSER_FOLDER})
ENDIF()


#
# Threading Building Blocks library
#
IF(NOT CMAKE_SYSTEM_NAME MATCHES "CYGWIN" AND
   NOT CMAKE_SYSTEM_NAME MATCHES "Windows")
  SET(FEATURE_THREADS_HAVE_BUNDLED TRUE)
  SET(TBB_FOLDER "${CMAKE_SOURCE_DIR}/bundled/tbb30_104oss")
  OPTION(DEAL_II_FORCE_BUNDLED_THREADS
    "Always use the bundled tbb library instead of an external one."
    OFF)
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


#
# UMFPACK, AMD and UFCONFIG:
#
SET(FEATURE_UMFPACK_HAVE_BUNDLED TRUE)
SET(UMFPACK_FOLDER "${CMAKE_SOURCE_DIR}/bundled/umfpack")
OPTION(DEAL_II_FORCE_BUNDLED_UMFPACK
  "Always use the bundled umfpack library instead of an external one."
  OFF)
MARK_AS_ADVANCED(
  atlas_LIBRARY
  blas_LIBRARY
  )

IF(FEATURE_UMFPACK_BUNDLED_CONFIGURED)
  ADD_SUBDIRECTORY(${UMFPACK_FOLDER}/UMFPACK/Source)
  ADD_SUBDIRECTORY(${UMFPACK_FOLDER}/AMD/Source)
ENDIF()

