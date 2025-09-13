## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2023 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------


##                                                                    ##
#            The cmake build system for the deal.II project            #
#                                                                      #
#   See doc/readme.html and doc/developers/cmake-internals.html for    #
#   further details on how to use the cmake build system of deal.II.   #
##                                                                    ##

########################################################################
#                                                                      #
#                            Configuration:                            #
#                                                                      #
########################################################################

#
# General configuration for cmake:
#
message(STATUS "This is CMake ${CMAKE_VERSION}")
message(STATUS "")

cmake_minimum_required(VERSION 3.13.4)
#
# We support all policy changes up to version 3.11.
#
cmake_policy(VERSION 3.13.4)

if(POLICY CMP0074)
  # Introduced in CMake 3.12: find_package(<PackageName>) also searches
  # <PackageName>_ROOT CMake and environment variables when set to NEW.
  cmake_policy(SET CMP0074 NEW)
endif()

if(POLICY CMP0075)
  # Introduced in CMake 3.12: Use CMAKE_REQUIRED_LIBRARIES also in include
  # file checks. We set the policy to NEW explicitly in order to avoid
  # spurious configure warnings.
  cmake_policy(SET CMP0075 NEW)
endif()

if(POLICY CMP0167)
  # Introduced in CMake 3.30: When set to true then CMake uses the Boost
  # CMake Config installed with Boost instead of using the FindBoost.cmake
  # module installed with CMake.
  cmake_policy(SET CMP0167 NEW)
endif()

if(POLICY CMP0177)
  # Introduced in CMAKE 3:31: Normalize DESTINATION paths so that they do
  # not contain "." or "..".
  cmake_policy(SET CMP0177 NEW)
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)

#
# Load all macros:
#
file(GLOB _macro_files "cmake/macros/*.cmake")
message(STATUS "Include ${CMAKE_SOURCE_DIR}/cmake/setup_external_macros.cmake")
include(${CMAKE_SOURCE_DIR}/cmake/setup_external_macros.cmake)
foreach(_file ${_macro_files})
  message(STATUS "Include ${_file}")
  include(${_file})
endforeach()

#
# Avoid verbose "Up-to-date" status information during installation:
#
set_if_empty(CMAKE_INSTALL_MESSAGE "LAZY")

#
# Check for the existence of various optional folders:
#
if(EXISTS ${CMAKE_SOURCE_DIR}/bundled/CMakeLists.txt)
  set(DEAL_II_HAVE_BUNDLED_DIRECTORY TRUE)
endif()

if(EXISTS ${CMAKE_SOURCE_DIR}/doc/CMakeLists.txt)
  set(DEAL_II_HAVE_DOC_DIRECTORY TRUE)
endif()

if(EXISTS ${CMAKE_SOURCE_DIR}/tests/CMakeLists.txt)
  set(DEAL_II_HAVE_TESTS_DIRECTORY TRUE)
endif()

#
# We have to initialize some cached variables before PROJECT is called, so
# do it at this point:
#
verbose_include(${CMAKE_SOURCE_DIR}/cmake/setup_cached_variables.cmake)

#
# Now, set the project and set up the rest:
#
project(deal.II CXX C)
enable_language_optional(Fortran)

verbose_include(${CMAKE_SOURCE_DIR}/cmake/setup_deal_ii.cmake)

verbose_include(${CMAKE_SOURCE_DIR}/cmake/setup_compiler_flags.cmake)

#
# Include information about bundled libraries:
#
if(DEAL_II_HAVE_BUNDLED_DIRECTORY)
  verbose_include(${CMAKE_SOURCE_DIR}/bundled/setup_bundled.cmake)
endif()

#
# Run all system checks:
#
file(GLOB _check_files "cmake/checks/*.cmake")
list(SORT _check_files)
foreach(_file ${_check_files})
  verbose_include(${_file})
endforeach()

# include setup_finalize.cmake before any target definitions
verbose_include(${CMAKE_SOURCE_DIR}/cmake/setup_finalize.cmake)

#
# Feature configuration:
#
file(GLOB _configure_files "cmake/configure/configure_*.cmake")
list(SORT _configure_files) # make sure to include in alphabetical order
foreach(_file ${_configure_files})
  verbose_include(${_file})
endforeach()

verbose_include(${CMAKE_SOURCE_DIR}/cmake/setup_sanity_checks.cmake)
verbose_include(${CMAKE_SOURCE_DIR}/cmake/setup_cpack.cmake)
verbose_include(${CMAKE_SOURCE_DIR}/cmake/setup_custom_targets.cmake)

########################################################################
#                                                                      #
#                     Compilation and installation:                    #
#                                                                      #
########################################################################

message(STATUS "")
message(STATUS "Configuring done. Proceed to target definitions now.")

set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

add_subdirectory(cmake/scripts)
add_subdirectory(include)

if(DEAL_II_HAVE_DOC_DIRECTORY)
  add_subdirectory(doc) # has to be included after include
endif()

if(DEAL_II_HAVE_BUNDLED_DIRECTORY)
  add_subdirectory(bundled)
endif()

# Now that we have collected information about the configuration
# and the bundled directory, set up the compilation of the
# library from the files in the source/ directory, as well as the
# C++20 module:
add_subdirectory(source)
add_subdirectory(module)

# Then also set up installation and other targets.
add_subdirectory(cmake/config)

# Now also configure installation of the examples. This relies on
# the correctly set up library, so needs to come after configuration
# in source/ and cmake/config/
add_subdirectory(examples)

add_subdirectory(contrib/utilities)

if(DEAL_II_HAVE_TESTS_DIRECTORY)
  add_subdirectory(tests)
endif()

add_subdirectory(contrib/python-bindings) # has to be included after tests

verbose_include(${CMAKE_SOURCE_DIR}/cmake/setup_write_config.cmake)

#
# And finally, print the configuration:
#
file(READ ${CMAKE_BINARY_DIR}/summary.log DEAL_II_LOG_SUMMARY)
message("${DEAL_II_LOG_SUMMARY}")
