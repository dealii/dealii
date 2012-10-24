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


##                                                                       ##
#                                                                         #
#             The cmake build system for the deal.II project              #
#                                                                         #
#   See README.cmake for further details on how to use the cmake build    #
#   system of deal.II.                                                    #
#                                                                         #
##                                                                       ##


###########################################################################
#                                                                         #
#                             Configuration:                              #
#                                                                         #
###########################################################################

SET(VERSION "8.0.pre")

#
# General configuration for cmake:
#
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

SET(CMAKE_INCLUDE_CURRENT_DIR true)

LIST(APPEND CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/cmake/
  ${CMAKE_SOURCE_DIR}/cmake/modules/
  )

#
# Load all macros:
#
FILE(GLOB _macro_files "cmake/macros/*.cmake")
FOREACH(_file ${_macro_files})
  MESSAGE(STATUS "Include ${_file}")
  INCLUDE(${_file})
ENDFOREACH()

INCLUDE(setup_external_macros)

#
# Check whether the doc or bundled folder is available:
#
IF(EXISTS ${CMAKE_SOURCE_DIR}/bundled/CMakeLists.txt)
  SET(DEAL_II_WITH_BUNDLED_DIRECTORY TRUE)
ENDIF()

IF(EXISTS ${CMAKE_SOURCE_DIR}/doc/CMakeLists.txt)
  SET(DEAL_II_WITH_DOC_DIRECTORY TRUE)
ENDIF()

#
# We have to initialize some cached variables before PROJECT is called, so
# do it at this point:
#
INCLUDE(setup_cached_variables)

#
# Now, set the project and setup the rest:
#
PROJECT(deal.II)
INCLUDE(setup_deal_ii)

INCLUDE(setup_compiler_flags)

#
# Run all system checks:
#
FILE(GLOB _check_files "cmake/checks/*.cmake")
FOREACH(_file ${_check_files})
  MESSAGE(STATUS "Include ${_file}")
  INCLUDE(${_file})
ENDFOREACH()

#
# Include information about bundled libraries:
#
IF(DEAL_II_WITH_BUNDLED_DIRECTORY)
  INCLUDE(bundled/CMakeLists.txt)
ENDIF()

#
# Feature configuration:
#
FILE(GLOB _configure_files "cmake/configure/configure_*.cmake")
# Make sure the configure_*.cmake files are included in alphabetical order
# to fullfill some dependencies...
LIST(SORT _configure_files)
FOREACH(_file ${_configure_files})
  MESSAGE(STATUS "Include ${_file}")
  INCLUDE(${_file})
ENDFOREACH()

#
# Finalize the configuration:
#
INCLUDE(setup_finalize)


###########################################################################
#                                                                         #
#                      Compilation and installation:                      #
#                                                                         #
###########################################################################

#
# Build the toolchain and prepare the includes:
#
INCLUDE(config/setup_template-arguments.in.cmake)
ADD_SUBDIRECTORY(scripts)
ADD_SUBDIRECTORY(include)

#
# Compile and install the library:
#
IF(DEAL_II_WITH_BUNDLED_DIRECTORY)
  ADD_SUBDIRECTORY(bundled)
ENDIF()
ADD_SUBDIRECTORY(source)

#
# Prepare the remaining directories for installation
#
ADD_SUBDIRECTORY(config)
IF(DEAL_II_WITH_DOC_DIRECTORY)
  ADD_SUBDIRECTORY(doc)
ENDIF()
ADD_SUBDIRECTORY(examples)
ADD_SUBDIRECTORY(contrib)

#
# Define some custom targets for convenience
#
INCLUDE(setup_custom_targets)

#
# And finally, print the configuration:
#
FILE(READ ${CMAKE_BINARY_DIR}/summary.log DEAL_II_LOG_SUMMARY)
MESSAGE("${DEAL_II_LOG_SUMMARY}")

