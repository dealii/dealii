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

#
# Set up deal.II specific definitions
#
# This file defines a long list of uncached variables, used throughout the
# configuration to determine paths, locations and names. Some linkage and
# crosscompilation setup happens also in here.
#
# Definitions marked with *) can be overridden by defining them to cache
# prior to the call of this file. This is done with the help of the
# SET_IF_EMPTY macro.
#
# General information about deal.II:
#
#     DEAL_II_PACKAGE_NAME            *)
#     DEAL_II_PACKAGE_VERSION         *)
#     DEAL_II_PACKAGE_VENDOR          *)
#     DEAL_II_PACKAGE_DESCRIPTION     *)
#     DEAL_II_VERSION_MAJOR
#     DEAL_II_VERSION_MINOR
#     DEAL_II_VERSION_SUBMINOR
#     DEAL_II_VERSION
#
# Information about paths, install locations and names:

#     DEAL_II_TARGET_NAME             *)
#
#     DEAL_II_PROJECT_CONFIG_NAME     *)
#     DEAL_II_BASE_NAME               *)
#     DEAL_II_DEBUG_SUFFIX            *)
#     DEAL_II_RELEASE_SUFFIX          *)
#
#     DEAL_II_EXECUTABLE_RELDIR       *)
#     DEAL_II_INCLUDE_RELDIR          *)
#     DEAL_II_LIBRARY_RELDIR          *)
#     DEAL_II_PKGCONF_RELDIR          *)
#     DEAL_II_PROJECT_CONFIG_RELDIR   *)
#     DEAL_II_SHARE_RELDIR            *)
#     DEAL_II_DOCREADME_RELDIR        *)
#     DEAL_II_DOCHTML_RELDIR          *)
#     DEAL_II_EXAMPLES_RELDIR         *)
#
#     DEAL_II_BUILD_TYPES
#     DEAL_II_LIST_SUFFIXES
#     DEAL_II_STRING_SUFFIXES
#
# *)  Can be overwritten by the command line via -D<...>
#

########################################################################
#                                                                      #
#                  General information about deal.II:                  #
#                                                                      #
########################################################################

set_if_empty(DEAL_II_PACKAGE_NAME "deal.II")

set_if_empty(DEAL_II_PACKAGE_VENDOR
  "The deal.II Authors <https://www.dealii.org/>"
  )
set_if_empty(DEAL_II_PACKAGE_DESCRIPTION
  "Library for solving partial differential equations with the finite element method"
  )

file(STRINGS "${CMAKE_SOURCE_DIR}/VERSION" _version LIMIT_COUNT 1)
set_if_empty(DEAL_II_PACKAGE_VERSION "${_version}")

#
# We expect a version number of the form "X.Y.Z" or "X.Y.Z-bla", where X, Y, Z
# are always numbers and bla is a short string ("pre", "rc0", "rc1", etc.).
#
string(REGEX REPLACE "^([0-9]+)\\..*" "\\1"
  DEAL_II_VERSION_MAJOR "${DEAL_II_PACKAGE_VERSION}"
  )
string(REGEX REPLACE "^[0-9]+\\.([0-9]+).*" "\\1"
  DEAL_II_VERSION_MINOR "${DEAL_II_PACKAGE_VERSION}"
  )
string(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1"
  DEAL_II_VERSION_SUBMINOR "${DEAL_II_PACKAGE_VERSION}"
  )

set(DEAL_II_VERSION ${DEAL_II_VERSION_MAJOR}.${DEAL_II_VERSION_MINOR}.${DEAL_II_VERSION_SUBMINOR})


########################################################################
#                                                                      #
#         Information about paths, install locations and names:        #
#                                                                      #
########################################################################

string(REPLACE "." "" _name "${DEAL_II_PACKAGE_NAME}")
string(TOLOWER "${_name}" _name)
set_if_empty(DEAL_II_TARGET_NAME "${_name}")

set_if_empty(DEAL_II_PROJECT_CONFIG_NAME "${DEAL_II_PACKAGE_NAME}")

string(REPLACE "." "_" _base_name "${DEAL_II_PACKAGE_NAME}")
set_if_empty(DEAL_II_BASE_NAME "${_base_name}")
set_if_empty(DEAL_II_DEBUG_SUFFIX ".g")
set_if_empty(DEAL_II_RELEASE_SUFFIX "")

#
# Try to obey the FSHS as close as possible ...
#
set_if_empty(DEAL_II_EXECUTABLE_RELDIR "bin")
set_if_empty(DEAL_II_INCLUDE_RELDIR "include")
set_if_empty(DEAL_II_LIBRARY_RELDIR "lib${LIB_SUFFIX}")
set_if_empty(DEAL_II_PKGCONF_RELDIR "${DEAL_II_LIBRARY_RELDIR}/pkgconfig")
set_if_empty(DEAL_II_PROJECT_CONFIG_RELDIR "${DEAL_II_LIBRARY_RELDIR}/cmake/${DEAL_II_PROJECT_CONFIG_NAME}")
set_if_empty(DEAL_II_SHARE_RELDIR "share/${DEAL_II_PACKAGE_NAME}")
#
# ... but install the documentation into prominent places:
#
set_if_empty(DEAL_II_DOCREADME_RELDIR "./")
set_if_empty(DEAL_II_DOCHTML_RELDIR "doc")
set_if_empty(DEAL_II_EXAMPLES_RELDIR "examples")

if(CMAKE_BUILD_TYPE MATCHES "Debug")
  list(APPEND DEAL_II_BUILD_TYPES "DEBUG")
endif()

if(CMAKE_BUILD_TYPE MATCHES "Release")
  list(APPEND DEAL_II_BUILD_TYPES "RELEASE")
endif()

set(DEAL_II_LIST_SUFFIXES
  LIBRARIES LIBRARIES_RELEASE LIBRARIES_DEBUG
  TARGETS TARGETS_RELEASE TARGETS_DEBUG
  INCLUDE_DIRS
  DEFINITIONS DEFINITIONS_RELEASE DEFINITIONS_DEBUG
  )

set(DEAL_II_STRING_SUFFIXES
  CXX_FLAGS CXX_FLAGS_RELEASE CXX_FLAGS_DEBUG
  LINKER_FLAGS LINKER_FLAGS_RELEASE LINKER_FLAGS_DEBUG
  EXECUTABLE
  )

#
# Disable platform introspection when cross compiling
#

if(CMAKE_CROSSCOMPILING)
  set(DEAL_II_ALLOW_PLATFORM_INTROSPECTION OFF CACHE BOOL "" FORCE)
endif()
