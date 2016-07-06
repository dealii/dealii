## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2016 by the deal.II authors
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


##                                                                    ##
#            The cmake build system for the deal.II project            #
#                                                                      #
#    See doc/readme.html and doc/development/cmake.html for further    #
#    details on how to use the cmake build system of deal.II.          #
##                                                                    ##

########################################################################
#                                                                      #
#                            Configuration:                            #
#                                                                      #
########################################################################

#
# General configuration for cmake:
#
MESSAGE(STATUS "This is CMake ${CMAKE_VERSION}")
MESSAGE(STATUS "")

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

#
# We support all policy changes up to version 3.1.0. Thus, explicitly set
# all policies CMP0001 - CMP0054 to new for version 3.1 (and later) to
# avoid some unnecessary warnings.
#
IF( "${CMAKE_VERSION}" VERSION_EQUAL "3.1" OR
    "${CMAKE_VERSION}" VERSION_GREATER "3.1" )
  CMAKE_POLICY(VERSION 3.1.0)
ENDIF()

IF(POLICY CMP0037)
  # allow to override "test" target for quick tests
  CMAKE_POLICY(SET CMP0037 OLD)
ENDIF()

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)

#
# Load all macros:
#
FILE(GLOB _macro_files "cmake/macros/*.cmake")
MESSAGE(STATUS "Include ${CMAKE_SOURCE_DIR}/cmake/setup_external_macros.cmake")
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/setup_external_macros.cmake)
FOREACH(_file ${_macro_files})
  MESSAGE(STATUS "Include ${_file}")
  INCLUDE(${_file})
ENDFOREACH()

#
# Check for the existence of various optional folders:
#
IF(EXISTS ${CMAKE_SOURCE_DIR}/bundled/CMakeLists.txt)
  SET(DEAL_II_HAVE_BUNDLED_DIRECTORY TRUE)
ENDIF()

IF(EXISTS ${CMAKE_SOURCE_DIR}/doc/CMakeLists.txt)
  SET(DEAL_II_HAVE_DOC_DIRECTORY TRUE)
ENDIF()

IF(EXISTS ${CMAKE_SOURCE_DIR}/tests/CMakeLists.txt)
  SET(DEAL_II_HAVE_TESTS_DIRECTORY TRUE)
ENDIF()

#
# We have to initialize some cached variables before PROJECT is called, so
# do it at this point:
#
VERBOSE_INCLUDE(${CMAKE_SOURCE_DIR}/cmake/setup_cached_variables.cmake)

#
# Now, set the project and set up the rest:
#
PROJECT(deal.II CXX C)
ENABLE_LANGUAGE_OPTIONAL(Fortran)

VERBOSE_INCLUDE(${CMAKE_SOURCE_DIR}/cmake/setup_deal_ii.cmake)

VERBOSE_INCLUDE(${CMAKE_SOURCE_DIR}/cmake/setup_compiler_flags.cmake)

#
# Include information about bundled libraries:
#
IF(DEAL_II_HAVE_BUNDLED_DIRECTORY)
  VERBOSE_INCLUDE(${CMAKE_SOURCE_DIR}/bundled/setup_bundled.cmake)
ENDIF()

#
# Run all system checks:
#
FILE(GLOB _check_files "cmake/checks/*.cmake")
LIST(SORT _check_files)
FOREACH(_file ${_check_files})
  VERBOSE_INCLUDE(${_file})
ENDFOREACH()

#
# Feature configuration:
#
FILE(GLOB _configure_files "cmake/configure/configure_*.cmake")
LIST(SORT _configure_files) # make sure to include in alphabetical order
FOREACH(_file ${_configure_files})
  VERBOSE_INCLUDE(${_file})
ENDFOREACH()

#
# Finalize the configuration:
#
VERBOSE_INCLUDE(${CMAKE_SOURCE_DIR}/cmake/setup_cpack.cmake)
VERBOSE_INCLUDE(${CMAKE_SOURCE_DIR}/cmake/setup_custom_targets.cmake)
VERBOSE_INCLUDE(${CMAKE_SOURCE_DIR}/cmake/setup_finalize.cmake)
VERBOSE_INCLUDE(${CMAKE_SOURCE_DIR}/cmake/setup_write_config.cmake)

########################################################################
#                                                                      #
#                     Compilation and installation:                    #
#                                                                      #
########################################################################

MESSAGE(STATUS "")
MESSAGE(STATUS "Configuring done. Proceed to target definitions now.")

ADD_SUBDIRECTORY(cmake/scripts)
ADD_SUBDIRECTORY(include)

IF(DEAL_II_HAVE_DOC_DIRECTORY)
  ADD_SUBDIRECTORY(doc) # has to be included after include
ENDIF()

IF(DEAL_II_HAVE_BUNDLED_DIRECTORY)
  ADD_SUBDIRECTORY(bundled)
ENDIF()
ADD_SUBDIRECTORY(source) # has to be included after bundled

ADD_SUBDIRECTORY(cmake/config) # has to be included after source
ADD_SUBDIRECTORY(examples)

IF(DEAL_II_HAVE_TESTS_DIRECTORY)
  ADD_SUBDIRECTORY(tests)
ENDIF()

#
# And finally, print the configuration:
#
FILE(READ ${CMAKE_BINARY_DIR}/summary.log DEAL_II_LOG_SUMMARY)
MESSAGE("${DEAL_II_LOG_SUMMARY}")
