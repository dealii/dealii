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
#             The cmake build system for the deal.II project
#

###########################################################################
##                                                                       ##
##  See README.cmake for further details on how to use the cmake build   ##
##  system of deal.II.                                                   ##
##                                                                       ##
###########################################################################

#
# The main purpose of this file is to set all the configuration options and
# then call a lot of modules residing under ./cmake/*/ that do the
# hard work.
#


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
FILE(GLOB macro_files "cmake/macros/*.cmake")
FOREACH(file ${macro_files})
  MESSAGE(STATUS "Include ${file}")
  INCLUDE(${file})
ENDFOREACH()

INCLUDE(setup_external_macros)

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
FILE(GLOB check_files "cmake/checks/*.cmake")
FOREACH(file ${check_files})
  MESSAGE(STATUS "Include ${file}")
  INCLUDE(${file})
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
FILE(GLOB configure_files "cmake/configure/configure_*.cmake")
# Make sure the configure_*.cmake files are included in alphabetical order
# to fullfill some dependencies...
LIST(SORT configure_files)
FOREACH(file ${configure_files})
  MESSAGE(STATUS "Include ${file}")
  INCLUDE(${file})
ENDFOREACH()

#
# Finalize the configuration:
#
INCLUDE(setup_finalize)


###########################################################################
#                                                                         #
#                    Compile and install the library:                     #
#                                                                         #
###########################################################################

#
# Build the toolchain and prepare the includes:
#
ADD_SUBDIRECTORY(config)
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
# Prepare some directories for installation
#
IF(DEAL_II_WITH_DOC_DIRECTORY)
  ADD_SUBDIRECTORY(doc)
ENDIF()
ADD_SUBDIRECTORY(examples)
ADD_SUBDIRECTORY(contrib)
IF(DEAL_II_COMPONENT_COMPAT_FILES)
  INCLUDE(config/setup_make.global_options.cmake)
ENDIF()


#
# Define some custom targets for convenience
#
INCLUDE(setup_custom_targets)

