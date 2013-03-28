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
# Set up deal.II specific definitions
#
# This file defines a long list of uncached variables, used throughout the
# configuration to determine paths, locations and names.
#
# Definitions marked with *) can be overriden by defining them to cache
# prior to the call of this file. This is done with the help of the
# SET_IF_EMPTY macro.
#
# General information about deal.II:
#
#     DEAL_II_PACKAGE_NAME            *)
#     DEAL_II_PACKAGE_VERSION         *)
#     DEAL_II_VERSION_MAJOR
#     DEAL_II_VERSION_MINOR
#     DEAL_II_VERSION
#
# Information about paths, install locations and names:
#
#     DEAL_II_PROJECT_CONFIG_NAME     *)
#     DEAL_II_BASE_NAME               *)
#     DEAL_II_DEBUG_SUFFIX            *)
#     DEAL_II_RELEASE_SUFFIX          *)
#
#     DEAL_II_PATH                    *)
#     DEAL_II_CMAKE_MACROS_RELDIR     *)
#     DEAL_II_DOCUMENTATION_RELDIR    *)
#     DEAL_II_EXAMPLES_RELDIR         *)
#     DEAL_II_EXECUTABLE_RELDIR       *)
#     DEAL_II_INCLUDE_RELDIR          *)
#     DEAL_II_LIBRARY_RELDIR          *)
#     DEAL_II_PROJECT_CONFIG_RELDIR   *)
#
#     DEAL_II_INCLUDE_DIRS
#     DEAL_II_LIBRARY_DIRS
#
#     DEAL_II_BUILD_TYPES
#     DEAL_II_WITH_BUNDLED_DIRECTORY
#     DEAL_II_WITH_DOC_DIRECTORY
#
# *)  Can be overwritten by the command line via -D<...>
#

###########################################################################
#                                                                         #
#                   General information about deal.II:                    #
#                                                                         #
###########################################################################

SET_IF_EMPTY(DEAL_II_PACKAGE_NAME "deal.II")

SET_IF_EMPTY(DEAL_II_PACKAGE_VERSION "8.0.pre") # TODO: Get this value from somewhere else

STRING(REGEX REPLACE
  "^([0-9]+)\\..*" "\\1" DEAL_II_VERSION_MAJOR "${DEAL_II_PACKAGE_VERSION}"
  )
STRING(REGEX REPLACE
  "^[0-9]+\\.([0-9]+).*" "\\1" DEAL_II_VERSION_MINOR "${DEAL_II_PACKAGE_VERSION}"
  )
SET(DEAL_II_VERSION ${DEAL_II_VERSION_MAJOR}.${DEAL_II_VERSION_MINOR})


###########################################################################
#                                                                         #
#          Information about paths, install locations and names:          #
#                                                                         #
###########################################################################

SET(DEAL_II_PROJECT_CONFIG_NAME "${DEAL_II_PACKAGE_NAME}")

STRING(REPLACE "." "_" _base_name "${DEAL_II_PACKAGE_NAME}")
SET_IF_EMPTY(DEAL_II_BASE_NAME "${_base_name}")
SET_IF_EMPTY(DEAL_II_DEBUG_SUFFIX ".g")
SET_IF_EMPTY(DEAL_II_RELEASE_SUFFIX "")

SET(DEAL_II_PATH ${CMAKE_INSTALL_PREFIX})

IF(DEAL_II_COMPONENT_COMPAT_FILES)
  #
  # The good, old directory structure:
  #
  SET_IF_EMPTY(DEAL_II_CMAKE_MACROS_RELDIR "cmake/macros")
  SET_IF_EMPTY(DEAL_II_COMMON_RELDIR "common")
  SET_IF_EMPTY(DEAL_II_DOCUMENTATION_RELDIR "doc")
  SET_IF_EMPTY(DEAL_II_EXAMPLES_RELDIR "examples")
  SET_IF_EMPTY(DEAL_II_EXECUTABLE_RELDIR "bin")
  SET_IF_EMPTY(DEAL_II_INCLUDE_RELDIR "include")
  SET_IF_EMPTY(DEAL_II_LIBRARY_RELDIR "lib")
  SET_IF_EMPTY(DEAL_II_PROJECT_CONFIG_RELDIR "${DEAL_II_LIBRARY_RELDIR}/cmake/${DEAL_II_PROJECT_CONFIG_NAME}")
ELSE()
  #
  # IF DEAL_II_COMPONENT_COMPAT_FILES=off, we assume that we have to
  # obey the FSHS...
  #
  SET_IF_EMPTY(DEAL_II_CMAKE_MACROS_RELDIR "share/${DEAL_II_PACKAGE_NAME}/cmake/Macros")
  SET_IF_EMPTY(DEAL_II_COMMON_RELDIR "share/${DEAL_II_PACKAGE_NAME}/common")
  SET_IF_EMPTY(DEAL_II_DOCUMENTATION_RELDIR "share/doc/${DEAL_II_PACKAGE_NAME}/html")
  SET_IF_EMPTY(DEAL_II_EXAMPLES_RELDIR "share/doc/${DEAL_II_PACKAGE_NAME}/examples")
  SET_IF_EMPTY(DEAL_II_EXECUTABLE_RELDIR "bin")
  SET_IF_EMPTY(DEAL_II_INCLUDE_RELDIR "include")
  SET_IF_EMPTY(DEAL_II_LIBRARY_RELDIR "lib${LIB_SUFFIX}")
  SET_IF_EMPTY(DEAL_II_PROJECT_CONFIG_RELDIR "${DEAL_II_LIBRARY_RELDIR}/cmake/${DEAL_II_PROJECT_CONFIG_NAME}")
ENDIF()

LIST(APPEND DEAL_II_INCLUDE_DIRS
  "${CMAKE_INSTALL_PREFIX}/${DEAL_II_INCLUDE_RELDIR}"
  "${CMAKE_INSTALL_PREFIX}/${DEAL_II_INCLUDE_RELDIR}/deal.II"
  "${CMAKE_INSTALL_PREFIX}/${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled"
  )

LIST(APPEND DEAL_II_LIBRARY_DIRS
  "${CMAKE_INSTALL_PREFIX}/${DEAL_II_LIBRARY_RELDIR}"
  )

IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  LIST(APPEND DEAL_II_BUILD_TYPES "DEBUG")
ENDIF()

IF(CMAKE_BUILD_TYPE MATCHES "Release")
  LIST(APPEND DEAL_II_BUILD_TYPES "RELEASE")
ENDIF()


#
# Cleanup some files used for storing the names of all object targets that
# will be bundled to the deal.II library.
# (Right now, i.e. cmake 2.8.8, this is the only reliable way to get
# information into a global scope...)
#
FOREACH(_build ${DEAL_II_BUILD_TYPES})
  STRING(TOLOWER "${_build}" _build_lowercase)
  FILE(REMOVE
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/deal_ii_objects_${_build_lowercase}
    )
ENDFOREACH()

