#####
##
## Copyright (C) 2012, 2013 by the deal.II authors
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
#     DEAL_II_ALLOW_PLATFORM_INTROSPECTION
#     DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS
#     CMAKE_BUILD_TYPE
#     BUILD_SHARED_LIBS
#     CMAKE_INSTALL_RPATH_USE_LINK_PATH
#     CMAKE_C_FLAGS
#     CMAKE_CXX_FLAGS                   *)
#     CMAKE_SHARED_LINKER_FLAGS         *)
#     DEAL_II_CXX_FLAGS_DEBUG
#     DEAL_II_SHARED_LINKER_FLAGS_DEBUG
#     DEAL_II_CXX_FLAGS_RELEASE
#     DEAL_II_SHARED_LINKER_FLAGS_RELEASE
#
# *)  May also be set via environment variable (CFLAGS, CXXFLAGS, LDFLAGS)
#     (nonempty cached variable has precedence will not be overwritten by
#     environment)
#


###########################################################################
#                                                                         #
#                     General configuration options:                      #
#                                                                         #
###########################################################################

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
  OFF
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

IF("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  SET(CMAKE_INSTALL_PREFIX
    "${CMAKE_BINARY_DIR}"
    CACHE STRING
    "Install path prefix, prepended onto install directories."
    )
ENDIF()


###########################################################################
#                                                                         #
#                        Compilation and linking:                         #
#                                                                         #
###########################################################################

OPTION(DEAL_II_ALLOW_PLATFORM_INTROSPECTION
  "Allow platform introspection for CPU command set, SSE and AVX"
  ON
  )
MARK_AS_ADVANCED(DEAL_II_ALLOW_PLATFORM_INTROSPECTION)

OPTION(DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS
  "Configure sensible default CFLAGS and CXXFLAGS depending on platform, compiler and build target."
  ON
  )
MARK_AS_ADVANCED(DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS)

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

SET(BUILD_SHARED_LIBS "ON" CACHE BOOL
  "Build a shared library"
  )
MARK_AS_ADVANCED(BUILD_SHARED_LIBS)

#
# Set CMAKE_INSTALL_RPATH_USE_LINK_PATH to default to ON and promote to
# cache so that the user can see the value.
#
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH "ON" CACHE BOOL
  "Set the rpath of the library to the external link pathes on installation"
  )
MARK_AS_ADVANCED(CMAKE_INSTALL_RPATH_USE_LINK_PATH)

#
# Tell the user very prominently, that we're doing things differently w.r.t
# CMAKE_(C|CXX)_FLAGS_(DEBUG|RELEASE)
#
FOREACH(_flag
  CXX_FLAGS_RELEASE
  CXX_FLAGS_DEBUG
  )
  IF(NOT "${CMAKE_${_flag}}" STREQUAL "")
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
  )
  # Go away...
  SET(${_flag} "" CACHE INTERNAL "" FORCE)
ENDFOREACH()

#
# Read in CFLAGS, CXXFLAGS and LDFLAGS from environment
#
SET_IF_EMPTY(CMAKE_C_FLAGS "$ENV{CFLAGS}")
SET_IF_EMPTY(CMAKE_CXX_FLAGS "$ENV{CXXFLAGS}")
SET_IF_EMPTY(CMAKE_SHARED_LINKER_FLAGS "$ENV{LDFLAGS}")

#
# Set cached compiler flags to an empty string:
#
SET(DEAL_II_USED_FLAGS
  CMAKE_CXX_FLAGS
  DEAL_II_CXX_FLAGS_DEBUG
  DEAL_II_CXX_FLAGS_RELEASE
  CMAKE_SHARED_LINKER_FLAGS
  DEAL_II_SHARED_LINKER_FLAGS_DEBUG
  DEAL_II_SHARED_LINKER_FLAGS_RELEASE
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


###########################################################################
#                                                                         #
#                          Miscellanious setup:                           #
#                                                                         #
###########################################################################

GET_CMAKE_PROPERTY(_res VARIABLES)
FOREACH(_var ${_res})
  #
  # Rename WITH_* by DEAL_II_WITH_*
  #
  IF(_var MATCHES "^WITH_")
    SET(DEAL_II_${_var} ${${_var}} CACHE BOOL "" FORCE)
    UNSET(${_var} CACHE)
  ENDIF()

  #
  # Same for components:
  #
  IF(_var MATCHES "^COMPONENT_")
    SET(DEAL_II_${_var} ${${_var}} CACHE BOOL "" FORCE)
    UNSET(${_var} CACHE)
  ENDIF()
  IF(_var MATCHES "^(COMPAT_FILES|DOCUMENTATION|EXAMPLES|MESH_CONVERTER|PARAMETER_GUI)")
    SET(DEAL_II_COMPONENT_${_var} ${${_var}} CACHE BOOL "" FORCE)
    UNSET(${_var} CACHE)
  ENDIF()

  #
  # If DEAL_II_FORCE_AUTODETECTION is set undefine all feature toggles
  # DEAL_II_WITH_* prior to configure:
  #
  IF(DEAL_II_FORCE_AUTODETECTION AND _var MATCHES "^DEAL_II_WITH_")
    UNSET(${_var} CACHE)
  ENDIF()
ENDFOREACH()
