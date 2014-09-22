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

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

IF(POLICY CMP0026)
  # enable target LOCATION property
  CMAKE_POLICY(SET CMP0026 OLD)
ENDIF()

IF(POLICY CMP0037)
  # allow to override "test" target for quick tests
  CMAKE_POLICY(SET CMP0037 OLD)
ENDIF()

MESSAGE(STATUS "This is CMake ${CMAKE_VERSION}")

LIST(APPEND CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/cmake/
  ${CMAKE_SOURCE_DIR}/cmake/modules/
  )

#
# Load all macros:
#
INCLUDE(setup_external_macros)
FILE(GLOB _macro_files "cmake/macros/*.cmake")
FOREACH(_file ${_macro_files})
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
INCLUDE(setup_cached_variables)

#
# Now, set the project and set up the rest:
#
PROJECT(deal.II CXX C)
ENABLE_LANGUAGE_OPTIONAL(Fortran)

INCLUDE(setup_deal_ii)

INCLUDE(setup_compiler_flags)

#
# Include information about bundled libraries:
#
IF(DEAL_II_HAVE_BUNDLED_DIRECTORY)
  SET(_file ${CMAKE_SOURCE_DIR}/bundled/setup_bundled.cmake)
  MESSAGE(STATUS "")
  MESSAGE(STATUS "Include ${_file}")
  INCLUDE(${_file})
ENDIF()

#
# Run all system checks:
#
FILE(GLOB _check_files "cmake/checks/*.cmake")
LIST(SORT _check_files)
FOREACH(_file ${_check_files})
  MESSAGE(STATUS "")
  MESSAGE(STATUS "Include ${_file}")
  INCLUDE(${_file})
ENDFOREACH()

#
# Feature configuration:
#
FILE(GLOB _configure_files "cmake/configure/configure_*.cmake")
LIST(SORT _configure_files) # make sure to include in alphabetical order
FOREACH(_file ${_configure_files})
  MESSAGE(STATUS "")
  MESSAGE(STATUS "Include ${_file}")
  INCLUDE(${_file})
ENDFOREACH()

#
# Finalize the configuration:
#
INCLUDE(setup_custom_targets)
INCLUDE(setup_finalize)
INCLUDE(setup_write_config)


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
ADD_SUBDIRECTORY(contrib) # has to be included after source
ADD_SUBDIRECTORY(examples)

IF(DEAL_II_HAVE_TESTS_DIRECTORY)
  ADD_SUBDIRECTORY(tests)
ENDIF()

#
# And finally, print the configuration:
#
FILE(READ ${CMAKE_BINARY_DIR}/summary.log DEAL_II_LOG_SUMMARY)
MESSAGE("${DEAL_II_LOG_SUMMARY}")
