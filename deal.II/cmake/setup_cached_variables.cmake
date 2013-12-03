## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2012 - 2013 by the deal.II authors
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

#
# Setup cached variables (prior to the PROJECT(deal.II) call)
#
# This file sets up the following cached Options:
#
# General configuration options:
#
#     DEAL_II_ALLOW_AUTODETECTION
#     DEAL_II_ALLOW_BUNDLED
#     DEAL_II_COMPONENT_COMPAT_FILES
#     DEAL_II_COMPONENT_DOCUMENTATION
#     DEAL_II_COMPONENT_EXAMPLES
#     DEAL_II_COMPONENT_MESH_CONVERTER
#     DEAL_II_COMPONENT_PARAMETER_GUI
#     DEAL_II_FORCE_AUTODETECTION
#
# Options regarding compilation and linking:
#
#     CMAKE_BUILD_TYPE
#     DEAL_II_ALLOW_PLATFORM_INTROSPECTION
#     DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS
#     DEAL_II_SETUP_COVERAGE
#     BUILD_SHARED_LIBS
#     DEAL_II_PREFER_STATIC_LIBS
#     DEAL_II_STATIC_EXECUTABLE
#     CMAKE_INSTALL_RPATH_USE_LINK_PATH
#     CMAKE_CXX_FLAGS                   *)
#     DEAL_II_LINKER_FLAGS              *)
#     DEAL_II_CXX_FLAGS_DEBUG
#     DEAL_II_LINKER_FLAGS_DEBUG
#     DEAL_II_CXX_FLAGS_RELEASE
#     DEAL_II_LINKER_FLAGS_RELEASE
#     DEAL_II_WITH_64BIT_INDICES
#
# *)  May also be set via environment variable (CXXFLAGS, LDFLAGS)
#     (a nonempty cached variable has precedence and will not be
#     overwritten by environment)
#


########################################################################
#                                                                      #
#                    General configuration options:                    #
#                                                                      #
########################################################################

If(DEAL_II_HAVE_BUNDLED_DIRECTORY)
  OPTION(DEAL_II_ALLOW_BUNDLED
    "Allow the use of libraries bundled with the source tarball. (DEAL_II_FORCE_BUNDLED* will overwrite this option.)"
    ON
    )
ENDIF()

OPTION(DEAL_II_COMPONENT_COMPAT_FILES
  "Enable installation of the example steps. This adds a COMPONENT \"compat_files\" to the build system."
  ON
  )

If(DEAL_II_HAVE_DOC_DIRECTORY)
  OPTION(DEAL_II_COMPONENT_DOCUMENTATION
    "Enable configuration, build and installation of the documentation. This adds a COMPONENT \"documentation\" to the build system."
    OFF
    )
ENDIF()

OPTION(DEAL_II_COMPONENT_EXAMPLES
  "Enable configuration and installation of the example steps. This adds a COMPONENT \"examples\" to the build system."
  ON
  )

OPTION(DEAL_II_COMPONENT_MESH_CONVERTER
  "Build and install the mesh_converter. This adds a COMPONENT \"mesh_converter\" to the build system."
  ON
  )

OPTION(DEAL_II_COMPONENT_PARAMETER_GUI
  "Build and install the parameter_gui. This adds a COMPONENT \"parameter_gui\" to the build system."
  OFF
  )

OPTION(DEAL_II_ALLOW_AUTODETECTION
  "Allow to automatically setup features by setting all undefined DEAL_II_WITH_* variables to ON or OFF"
  ON
  )

OPTION(DEAL_II_FORCE_AUTODETECTION
  "Force feature autodetection by undefining all DEAL_II_WITH_* variables prior to configure"
  OFF
  )


########################################################################
#                                                                      #
#                       Compilation and linking:                       #
#                                                                      #
########################################################################

#
# Setup CMAKE_BUILD_TYPE:
#
SET(CMAKE_BUILD_TYPE
  "DebugRelease"
  CACHE STRING
  "Choose the type of build, options are: Debug, Release and DebugRelease."
  )

#
# This is cruel, I know. But it is better to only have a known number of
# options for CMAKE_BUILD_TYPE...
#
IF( NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Release" AND
    NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" AND
    NOT "${CMAKE_BUILD_TYPE}" STREQUAL "DebugRelease" )
  MESSAGE(FATAL_ERROR
    "CMAKE_BUILD_TYPE does neither match Release, Debug, nor DebugRelease!"
    )
ENDIF()


OPTION(DEAL_II_ALLOW_PLATFORM_INTROSPECTION
  "Allow platform introspection for CPU command sets, SSE and AVX"
  ON
  )
MARK_AS_ADVANCED(DEAL_II_ALLOW_PLATFORM_INTROSPECTION)

OPTION(DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS
  "Configure sensible default CFLAGS and CXXFLAGS depending on platform, compiler and build target."
  ON
  )
MARK_AS_ADVANCED(DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS)

OPTION(DEAL_II_SETUP_COVERAGE
  "Setup debug compiler flags to provide additional test coverage information. Currently only gprof is supported."
  OFF
  )
MARK_AS_ADVANCED(DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS)

SET(BUILD_SHARED_LIBS "ON" CACHE BOOL
  "Build a shared library"
  )

OPTION(DEAL_II_PREFER_STATIC_LIBS
  "Prefer static libraries over dynamic libraries when searching for features and corresponding link interface"
  OFF
  )
MARK_AS_ADVANCED(DEAL_II_PREFER_STATIC_LIBS)

OPTION(DEAL_II_STATIC_EXECUTABLE
  "Provide a link interface that is suitable for static linkage of executables. Enabling this option forces BUILD_SHARED_LIBS=OFF and DEAL_II_PREFER_STATIC_LIBS=ON"
  OFF
  )
MARK_AS_ADVANCED(DEAL_II_STATIC_EXECUTABLE)

IF(DEAL_II_STATIC_EXECUTABLE)
  SET(BUILD_SHARED_LIBS "OFF" CACHE BOOL
    "Build a shared library"
    FORCE
    )
  SET(DEAL_II_PREFER_STATIC_LIBS "ON" CACHE BOOL
    "Prefer static libraries over dynamic libraries when searching for features and corresponding link interface"
    FORCE
    )
ENDIF()


#
# Set CMAKE_INSTALL_RPATH_USE_LINK_PATH to default to ON and promote to
# cache so that the user can see the value.
#
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH "ON" CACHE BOOL
  "Set the rpath of the library to the external link pathes on installation"
  )
MARK_AS_ADVANCED(CMAKE_INSTALL_RPATH_USE_LINK_PATH)


#
# Define the variable that defines whether we should use 32- or 64-bit
# global DoF indices.
#
OPTION(DEAL_II_WITH_64BIT_INDICES
  "If set to ON, then use 64-bit data types to represent global degree of freedom indices. The default is to OFF. You only want to set this to ON if you will solve problems with more than 2^31 (approximately 2 billion) unknowns. If set to ON, you also need to ensure that both Trilinos and/or PETSc support 64-bit indices."
  OFF
  )


#
# Tell the user very prominently, that we're doing things differently w.r.t
# CMAKE_(C|CXX)_FLAGS_(DEBUG|RELEASE)
#
FOREACH(_flag
  CXX_FLAGS_RELEASE
  CXX_FLAGS_DEBUG
  SHARED_LINKER_FLAGS_RELEASE
  SHARED_LINKER_FLAGS_DEBUG
  )
  IF(NOT "${CMAKE_${_flag}}" STREQUAL "")
    UNSET(${CMAKE_${_flag}} CACHE)
    MESSAGE(FATAL_ERROR
      "\nThe deal.II cmake build system does not use CMAKE_${_flag}.\n"
      "Use DEAL_II_${_flag}, instead!\n\n"
      )
  ENDIF()
ENDFOREACH()


#
# Hide all unused compiler flag variables:
#
FOREACH(_flag
  CMAKE_CXX_FLAGS_RELEASE
  CMAKE_CXX_FLAGS_DEBUG
  CMAKE_CXX_FLAGS_MINSIZEREL
  CMAKE_CXX_FLAGS_RELWITHDEBINFO
  CMAKE_C_FLAGS_RELEASE
  CMAKE_C_FLAGS_DEBUG
  CMAKE_C_FLAGS_MINSIZEREL
  CMAKE_C_FLAGS_RELWITHDEBINFO
  CMAKE_Fortran_FLAGS_RELEASE
  CMAKE_Fortran_FLAGS_DEBUG
  CMAKE_Fortran_FLAGS_MINSIZEREL
  CMAKE_Fortran_FLAGS_RELWITHDEBINFO
  CMAKE_SHARED_LINKER_FLAGS
  CMAKE_SHARED_LINKER_FLAGS_DEBUG
  CMAKE_SHARED_LINKER_FLAGS_MINSIZEREL
  CMAKE_SHARED_LINKER_FLAGS_RELEASE
  CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO
  )
  # Go away...
  SET(${_flag} "" CACHE INTERNAL "" FORCE)
ENDFOREACH()


#
# Set cached compiler flags to an empty string:
#
SET(DEAL_II_USED_FLAGS
  CMAKE_CXX_FLAGS
  DEAL_II_CXX_FLAGS_DEBUG
  DEAL_II_CXX_FLAGS_RELEASE
  DEAL_II_LINKER_FLAGS
  DEAL_II_LINKER_FLAGS_DEBUG
  DEAL_II_LINKER_FLAGS_RELEASE
  )

FOREACH(_flag ${DEAL_II_USED_FLAGS})
  #
  # Promote to cache:
  #
  SET(${_flag} "${${_flag}}" CACHE STRING
   "The user supplied cache variable will be appended _at the end_ of the auto generated ${_flag} variable"
   )

  #
  # Save the initial (cached) variable at this point and clear it.
  # ${flags}_SAVED will be appended to ${flags} in
  # setup_finalize.cmake (called at the end of the
  # main CMakeLists.txt file).
  #
  SET(${_flag}_SAVED "${${_flag}}")
  SET(${_flag} "")

  #
  # Mark these flags as advanced.
  #
  MARK_AS_ADVANCED(${_flag})
ENDFOREACH()


#
# Finally, read in CXXFLAGS and LDFLAGS from environment and prepend them
# to the saved variables:
#
SET(CMAKE_CXX_FLAGS_SAVED "$ENV{CXXFLAGS} ${CMAKE_CXX_FLAGS_SAVED}")
SET(DEAL_II_LINKER_FLAGS_SAVED "$ENV{LDFLAGS} ${DEAL_II_LINKER_FLAGS_SAVED}")
UNSET(ENV{CXXFLAGS})
UNSET(ENV{LDFLAGS})

#
# Also respect DEAL_II_CXX_FLAGS - it just too easy to accientally write
# DEAL_II_CXX_FLAGS instead of CMAKE_CXX_FLAGS... (this was a poor design
# choice, I know...)
#
IF(NOT "${DEAL_II_CXX_FLAGS}" STREQUAL "")
  MESSAGE(STATUS
    "Appending \${DEAL_II_CXX_FLAGS} to saved \${CMAKE_CXX_FLAGS}"
    )
  SET(CMAKE_CXX_FLAGS_SAVED "${CMAKE_CXX_FLAGS_SAVED} ${DEAL_II_CXX_FLAGS}")
ENDIF()



########################################################################
#                                                                      #
#                         Miscellaneous setup:                         #
#                                                                      #
########################################################################

#
# We do not support installation into the binary directory any more ("too
# much pain, not enough profit"):
#

IF("${CMAKE_BINARY_DIR}" STREQUAL "${CMAKE_INSTALL_PREFIX}")
  MESSAGE(FATAL_ERROR "
Error CMAKE_INSTALL_PREFIX is equal to CMAKE_BINARY_DIR.
It is not possible to install into the build directory. Please set
CMAKE_INSTALL_PREFIX to a designated install directory different than
CMAKE_BINARY_DIR.
(Please note that you can use deal.II directly out of a build directory
if this is what you tried to do.)
"
    )
ENDIF()


GET_CMAKE_PROPERTY(_res VARIABLES)
FOREACH(_var ${_res})
  #
  # Rename (ALLOW|WITH|FORCE|COMPONENT)_* by DEAL_II_(ALLOW|WITH|FORCE|COMPONENT)_*
  #
  FOREACH(_match ALLOW_ WITH_ FORCE_ COMPONENT_)
    IF(_var MATCHES "^${_match}")
      SET(DEAL_II_${_var} ${${_var}} CACHE BOOL "" FORCE)
      UNSET(${_var} CACHE)
    ENDIF()
  ENDFOREACH()

  #
  # Same for components:
  #
  IF(_var MATCHES "^(COMPAT_FILES|DOCUMENTATION|EXAMPLES|MESH_CONVERTER|PARAMETER_GUI)")
    SET(DEAL_II_COMPONENT_${_var} ${${_var}} CACHE BOOL "" FORCE)
    UNSET(${_var} CACHE)
  ENDIF()

  #
  # If DEAL_II_FORCE_AUTODETECTION is set undefine all feature toggles
  # DEAL_II_WITH_* prior to configure:
  #
  IF(DEAL_II_FORCE_AUTODETECTION AND _var MATCHES "^DEAL_II_WITH_"
     # Exclude FEATURES that do not represent external libraries:
     AND NOT _var MATCHES "^DEAL_II_WITH_64BIT_INDICES" )
    UNSET(${_var} CACHE)
  ENDIF()
ENDFOREACH()
