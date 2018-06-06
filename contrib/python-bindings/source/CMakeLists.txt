## ---------------------------------------------------------------------
##
## Copyright (C) 2016 - 2018 by the deal.II authors
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
# Find Python:
#

INCLUDE(FindPythonInterp)
INCLUDE(FindPythonLibs)

IF(FEATURE_BOOST_BUNDLED_CONFIGURED)
  MESSAGE(FATAL_ERROR
    "DEAL_II_COMPONENT_PYTHON_BINDINGS has unmet configuration requirements: "
    "Python bindings require an external boost library, but deal.II was "
    "configured with bundled boost."
    )
ENDIF()

#
# Unset Boost_Found and run the low level FindBOOST.cmake module again to
# pick up libboost_python.so
#
LIST(REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)
SET(Boost_FOUND)
# Use the low level _FIND_PACKAGE function instead of our wrapper
_FIND_PACKAGE(Boost 1.54 COMPONENTS python)
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)

IF(NOT Boost_FOUND)
  MESSAGE(FATAL_ERROR
    "DEAL_II_COMPONENT_PYTHON_BINDINGS has unmet configuration requirements: "
    "The external boost library does not provide Boost.Python"
    )
ENDIF()

#
# FIXME: Once finalized, reconsider moving this definitions into
# cmake/setup_dealii.cmake
#
# General information about deal.II:
#
#     PYTHON_BINDINGS_BASE_NAME       *)
#     PYTHON_BINDINGS_DEBUG_NAME      *)
#     PYTHON_BINDINGS_RELEASE_NAME    *)
#
# Information about paths, install locations and names:
#
#     DEAL_II_PYTHON_RELDIR           *) **)
#
# *)  Can be overwritten by the command line via -D<...>
# **) We do a best effort guess on site-packages location...
#

SET_IF_EMPTY(PYTHON_BINDINGS_BASE_NAME "PyDealII")
SET_IF_EMPTY(PYTHON_BINDINGS_DEBUG_NAME "Debug")
SET_IF_EMPTY(PYTHON_BINDINGS_RELEASE_NAME "Release")

SET_IF_EMPTY(DEAL_II_PYTHON_RELDIR
  "${DEAL_II_LIBRARY_RELDIR}/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages/${PYTHON_BINDINGS_BASE_NAME}"
  )

INCLUDE_DIRECTORIES(
  ${CMAKE_BINARY_DIR}/include/
  ${CMAKE_SOURCE_DIR}/include/
  ${DEAL_II_BUNDLED_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/contrib/python-bindings/include/
  SYSTEM
  ${DEAL_II_INCLUDE_DIRS}
  )

#
# Build a variant of the wrapper library for all configured build types in
# DEAL_II_BUILD_TYPES (usually PyDealIIDebug and PyDealIIRelease)
#

INCLUDE_DIRECTORIES(SYSTEM ${PYTHON_INCLUDE_DIRS})

SET(_src
  wrappers.cc
  export_cell_accessor.cc
  export_point.cc
  export_triangulation.cc
  cell_accessor_wrapper.cc
  point_wrapper.cc
  triangulation_wrapper.cc
  )

FOREACH(_build ${DEAL_II_BUILD_TYPES})
  STRING(TOLOWER ${_build} _build_lowercase)

  PYTHON_ADD_MODULE(PyDealII_${_build_lowercase} ${_src})

  SET_TARGET_PROPERTIES(PyDealII_${_build_lowercase} PROPERTIES
    OUTPUT_NAME "${PYTHON_BINDINGS_${_build}_NAME}"
    LINK_FLAGS "${DEAL_II_LINKER_FLAGS} ${DEAL_II_LINKER_FLAGS_${_build}}"
    LINKER_LANGUAGE "CXX"
    COMPILE_DEFINITIONS "${DEAL_II_DEFINITIONS};${DEAL_II_DEFINITIONS_${_build}}"
    COMPILE_FLAGS "${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${_build}}"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${DEAL_II_PYTHON_RELDIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${DEAL_II_PYTHON_RELDIR}"
    )
  IF(CMAKE_SYSTEM_NAME MATCHES "Darwin")
    SET_TARGET_PROPERTIES(PyDealII_${_build_lowercase} PROPERTIES
      MACOSX_RPATH OFF
      BUILD_WITH_INSTALL_RPATH OFF
      INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/${DEAL_II_PYTHON_RELDIR}"
      )
  ENDIF()

  TARGET_LINK_LIBRARIES(PyDealII_${_build_lowercase}
    ${DEAL_II_BASE_NAME}${DEAL_II_${_build}_SUFFIX}
    ${Boost_LIBRARIES}
    ${PYTHON_LIBRARIES}
    )

  EXPORT(TARGETS PyDealII_${_build_lowercase}
    FILE ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}Targets.cmake
    APPEND
    )

  INSTALL(TARGETS  PyDealII_${_build_lowercase}
    COMPONENT python_bindings
    EXPORT ${DEAL_II_PROJECT_CONFIG_NAME}Targets
    LIBRARY DESTINATION ${DEAL_II_PYTHON_RELDIR}
    ARCHIVE DESTINATION ${DEAL_II_PYTHON_RELDIR}
    )

  ADD_DEPENDENCIES(python_bindings PyDealII_${_build_lowercase})
ENDFOREACH()

#
# Copy python sources to binary directory and set up installation:
#

SET(_python_sources
  __init__.py
  )
FILE(COPY ${_python_sources}
  DESTINATION ${CMAKE_BINARY_DIR}/${DEAL_II_PYTHON_RELDIR}
  )
INSTALL(FILES ${_python_sources}
  DESTINATION ${DEAL_II_PYTHON_RELDIR}
  COMPONENT python_bindings
  )

