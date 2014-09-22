# ---------------------------------------------------------------------
##
## Copyright (C) 2013, 2014 by the deal.II authors
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

SET(BOOST_FOLDER "${CMAKE_SOURCE_DIR}/bundled/boost-1.56.0")


#
# Threading Building Blocks library
#

IF( NOT CMAKE_SYSTEM_NAME MATCHES "CYGWIN"
    AND NOT CMAKE_SYSTEM_NAME MATCHES "Windows"
    AND NOT CMAKE_SYSTEM_PROCESSOR MATCHES "ia64" )
  #
  # Cygwin is unsupported by tbb, Windows due to the way we compile tbb...
  #
  SET(FEATURE_THREADS_HAVE_BUNDLED TRUE)

  OPTION(DEAL_II_FORCE_BUNDLED_THREADS
    "Always use the bundled tbb library instead of an external one."
    OFF)

  SET(TBB_FOLDER "${CMAKE_SOURCE_DIR}/bundled/tbb41_20130401oss")
ENDIF()

#
# UMFPACK, AMD and UFCONFIG:
#

SET(FEATURE_UMFPACK_HAVE_BUNDLED TRUE)

OPTION(DEAL_II_FORCE_BUNDLED_UMFPACK
  "Always use the bundled umfpack library instead of an external one."
  OFF)

SET(UMFPACK_FOLDER "${CMAKE_SOURCE_DIR}/bundled/umfpack")


#
# muparser
#
SET(FEATURE_MUPARSER_HAVE_BUNDLED TRUE)

OPTION(DEAL_II_FORCE_BUNDLED_MUPARSER
  "Always use the bundled functionparser library instead of an external one."
  OFF)

SET(MUPARSER_FOLDER "${CMAKE_SOURCE_DIR}/bundled/muparser_v2_2_3/")
