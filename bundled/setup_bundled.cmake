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
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

#
# Export information about bundled library locations and do the actual
# setup of compilation targets and installation here:
#


#
# Boost C++ libraries
#

set(FEATURE_BOOST_HAVE_BUNDLED TRUE)

option(DEAL_II_FORCE_BUNDLED_BOOST
  "Always use the bundled boost library instead of an external one."
  OFF)

set(BOOST_FOLDER "${CMAKE_SOURCE_DIR}/bundled/boost-1.70.0")

macro(feature_boost_configure_bundled)
  set(BOOST_VERSION "1.70.0")
  set(BOOST_VERSION_MAJOR "1")
  set(BOOST_VERSION_MINOR "70")
  set(BOOST_VERSION_SUBMINOR "0")

  #
  # Add rt to the link interface as well, boost/chrono needs it.
  #
  if(NOT CMAKE_SYSTEM_NAME MATCHES "Windows")
    find_system_library(rt_LIBRARY NAMES rt)
    mark_as_advanced(rt_LIBRARY)
    if(NOT rt_LIBRARY MATCHES "-NOTFOUND")
      list(APPEND DEAL_II_LIBRARIES ${rt_LIBRARY})
    endif()
  endif()

  if(CMAKE_SYSTEM_NAME MATCHES "Windows")
    #
    # Bundled boost tries to (dl)open itself as a dynamic library on
    # Windows. Disable this undesired behavior by exporting
    # BOOST_ALL_NO_LIB on Windows platforms (for bundled boost).
    #
    list(APPEND DEAL_II_DEFINITIONS "BOOST_ALL_NO_LIB")
  endif()

  enable_if_supported(DEAL_II_WARNING_FLAGS "-Wno-unused-local-typedefs")

  list(APPEND DEAL_II_BUNDLED_INCLUDE_DIRS ${BOOST_FOLDER}/include)
endmacro()


#
# Kokkos
#

set(FEATURE_KOKKOS_HAVE_BUNDLED TRUE)

option(DEAL_II_FORCE_BUNDLED_KOKKOS
  "Always use the bundled Kokkos library instead of an external one."
  OFF)

set(KOKKOS_FOLDER "${CMAKE_SOURCE_DIR}/bundled/kokkos-3.7.00")

macro(feature_kokkos_configure_bundled)
  set(KOKKOS_VERSION "3.7.0")
  set(Kokkos_DEVICES "Serial")
  set(Kokkos_ARCH " ")

  list(APPEND DEAL_II_BUNDLED_INCLUDE_DIRS
    ${KOKKOS_FOLDER}/algorithms/src
    ${KOKKOS_FOLDER}/containers/src
    ${KOKKOS_FOLDER}/core/src
    ${KOKKOS_FOLDER}/simd/src
    ${KOKKOS_FOLDER}/tpls/desul/include
    )
endmacro()


#
# Taskflow
#

set(FEATURE_TASKFLOW_HAVE_BUNDLED TRUE)

option(DEAL_II_FORCE_BUNDLED_TASKFLOW
  "Always use the bundled taskflow header library instead of an external one."
  OFF)

set(TASKFLOW_FOLDER "${CMAKE_SOURCE_DIR}/bundled/taskflow-2.5.0")

macro(feature_taskflow_configure_bundled)
  set(TASKFLOW_VERSION "2.5.0")

  list(APPEND DEAL_II_BUNDLED_INCLUDE_DIRS ${TASKFLOW_FOLDER}/include)
endmacro()


#
# Threading Building Blocks library
#

if( NOT CMAKE_SYSTEM_NAME MATCHES "CYGWIN"
    AND NOT CMAKE_SYSTEM_NAME MATCHES "Windows" )
  #
  # Cygwin is unsupported by tbb, Windows due to the way we compile tbb...
  #
  set(FEATURE_TBB_HAVE_BUNDLED TRUE)

  option(DEAL_II_FORCE_BUNDLED_TBB
    "Always use the bundled tbb library instead of an external one."
    OFF)

  set(TBB_FOLDER "${CMAKE_SOURCE_DIR}/bundled/tbb-2018_U2")
endif()

macro(feature_tbb_configure_bundled)
  set(TBB_VERSION "2018.0")
  set(TBB_VERSION_MAJOR "2018")
  set(TBB_VERSION_MINOR "0")

  #
  # We have to disable a bunch of warnings:
  #
  enable_if_supported(DEAL_II_WARNING_FLAGS "-Wno-parentheses")

  #
  # tbb uses dlopen/dlclose, so link against libdl.so as well:
  #
  list(APPEND DEAL_II_LIBRARIES ${CMAKE_DL_LIBS})

  list(APPEND DEAL_II_BUNDLED_INCLUDE_DIRS ${TBB_FOLDER}/include)

  set(DEAL_II_TBB_WITH_ONEAPI FALSE)
endmacro()


#
# UMFPACK, AMD and UFCONFIG:
#

set(FEATURE_UMFPACK_HAVE_BUNDLED TRUE)

option(DEAL_II_FORCE_BUNDLED_UMFPACK
  "Always use the bundled umfpack library instead of an external one."
  OFF)

set(UMFPACK_FOLDER "${CMAKE_SOURCE_DIR}/bundled/umfpack")

macro(feature_umfpack_configure_bundled)
  set(UMFPACK_VERSION "5.0.2")
  set(UMFPACK_VERSION_MAJOR "5")
  set(UMFPACK_VERSION_MINOR "0")
  set(UMFPACK_VERSION_SUBMINOR "2")

  list(APPEND DEAL_II_BUNDLED_INCLUDE_DIRS
    ${UMFPACK_FOLDER}/UMFPACK/Include ${UMFPACK_FOLDER}/AMD/Include
    )
endmacro()


#
# muparser
#
set(FEATURE_MUPARSER_HAVE_BUNDLED TRUE)

option(DEAL_II_FORCE_BUNDLED_MUPARSER
  "Always use the bundled functionparser library instead of an external one."
  OFF)

set(MUPARSER_FOLDER "${CMAKE_SOURCE_DIR}/bundled/muparser_v2_3_3/")

macro(feature_muparser_configure_bundled)
  set(MUPARSER_VERSION "2.3.3")
  set(MUPARSER_VERSION_MAJOR "2")
  set(MUPARSER_VERSION_MINOR "3")
  set(MUPARSER_VERSION_SUBMINOR "3")

  list(APPEND DEAL_II_BUNDLED_INCLUDE_DIRS ${MUPARSER_FOLDER}/include)
endmacro()
