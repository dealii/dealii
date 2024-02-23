## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2016 - 2023 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

if(DEAL_II_COMPONENT_PYTHON_BINDINGS)

  message(STATUS "")
  message(STATUS "Setting up python bindings")

  #
  # Find Python:
  #
  find_package(Python3)
  set(PYTHON_VERSION_MAJOR ${Python3_VERSION_MAJOR})
  set(PYTHON_VERSION_MINOR ${Python3_VERSION_MINOR})
  set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE})
  include(FindPythonLibs)

  if(DEAL_II_FEATURE_BOOST_BUNDLED_CONFIGURED)
    message(FATAL_ERROR
      "DEAL_II_COMPONENT_PYTHON_BINDINGS has unmet configuration requirements: "
      "Python bindings require an external boost library, but deal.II was "
      "configured with bundled boost."
      )
  endif()

  #
  # As of 1.67, boost requires specifying the suffix for the python
  # component manually.
  #
  unset(Boost_FOUND)
  # On case-insensitive file systems, FindBOOST.cmake and FindBoost.cmake are indistinguishable.
  # Make sure FindBoost.cmake from CMake installation is found.
  list(REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)
  set(Boost_USE_DEBUG_RUNTIME FALSE)
  if(${BOOST_VERSION} VERSION_LESS 1.67)
    find_package(Boost 1.59 COMPONENTS python REQUIRED)
  else()
    find_package(Boost 1.67 COMPONENTS python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} REQUIRED)
  endif()
  list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)

  if(NOT Boost_FOUND)
    message(FATAL_ERROR
      "DEAL_II_COMPONENT_PYTHON_BINDINGS has unmet configuration requirements: "
      "The external boost library does not provide Boost.Python"
      )
  endif()

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

  set_if_empty(PYTHON_BINDINGS_BASE_NAME "PyDealII")
  set_if_empty(PYTHON_BINDINGS_DEBUG_NAME "Debug")
  set_if_empty(PYTHON_BINDINGS_RELEASE_NAME "Release")

  set_if_empty(DEAL_II_PYTHON_RELDIR
    "${DEAL_II_LIBRARY_RELDIR}/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages/${PYTHON_BINDINGS_BASE_NAME}"
    )

  add_subdirectory(source)

  add_subdirectory(tests)

  message(STATUS "Setting up python bindings - Done")
  message(STATUS "")

endif()
