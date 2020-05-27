## ---------------------------------------------------------------------
##
## Copyright (C) 2016 - 2020 by the deal.II authors
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

IF(DEAL_II_COMPONENT_PYTHON_BINDINGS)

  MESSAGE(STATUS "")
  MESSAGE(STATUS "Setting up python bindings")

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
  # As of 1.67, boost requires specifying the suffix for the python
  # component manually.
  #
  UNSET(Boost_FOUND)
  # On case-insensitive file systems, FindBOOST.cmake and FindBoost.cmake are indistinguishable.
  # Make sure FindBoost.cmake from CMake installation is found.
  LIST(REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)
  IF(${BOOST_VERSION} VERSION_LESS 1.67)
    _FIND_PACKAGE(Boost 1.59 COMPONENTS python REQUIRED)
  ELSE()
    _FIND_PACKAGE(Boost 1.67 COMPONENTS python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} REQUIRED)
  ENDIF()
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

  ADD_SUBDIRECTORY(source)

  ADD_SUBDIRECTORY(tests)

  MESSAGE(STATUS "Setting up python bindings - Done")
  MESSAGE(STATUS "")

ENDIF(DEAL_II_COMPONENT_PYTHON_BINDINGS)
