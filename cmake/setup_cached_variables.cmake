## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2017 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

#
# Set up cached variables (prior to the PROJECT(deal.II) call)
#
# This file sets up the following cached Options:
#
# General configuration options:
#
#     DEAL_II_ALLOW_AUTODETECTION
#     DEAL_II_ALLOW_BUNDLED
#     DEAL_II_COMPONENT_DOCUMENTATION
#     DEAL_II_COMPONENT_EXAMPLES
#     DEAL_II_COMPONENT_PACKAGE
#     DEAL_II_COMPONENT_PYTHON_BINDINGS
#     DEAL_II_FORCE_AUTODETECTION
#
# Options regarding compilation and linking:
#
#     CMAKE_BUILD_TYPE
#     DEAL_II_ALLOW_PLATFORM_INTROSPECTION
#     DEAL_II_WITH_SIMD
#     DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS
#     DEAL_II_SETUP_COVERAGE
#     DEAL_II_UNITY_BUILD
#     BUILD_SHARED_LIBS
#     DEAL_II_PREFER_STATIC_LIBS
#     DEAL_II_STATIC_EXECUTABLE
#     CMAKE_INSTALL_RPATH_USE_LINK_PATH
#     DEAL_II_CXX_FLAGS                    *)
#     DEAL_II_CXX_FLAGS_DEBUG
#     DEAL_II_CXX_FLAGS_RELEASE
#     DEAL_II_LINKER_FLAGS                 *)
#     DEAL_II_LINKER_FLAGS_DEBUG
#     DEAL_II_LINKER_FLAGS_RELEASE
#
# Components and miscellaneous options:
#
#     DEAL_II_WITH_64BIT_INDICES
#     DEAL_II_WITH_COMPLEX_VALUES
#     DEAL_II_DOXYGEN_USE_MATHJAX
#     DEAL_II_COMPILE_EXAMPLES
#     DEAL_II_CPACK_BUNDLE_NAME
#     DEAL_II_CPACK_EXTERNAL_LIBS
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

If(DEAL_II_HAVE_DOC_DIRECTORY)
  OPTION(DEAL_II_COMPONENT_DOCUMENTATION
    "Enable configuration, build and installation of the documentation. This adds a COMPONENT \"documentation\" to the build system."
    OFF
    )
  LIST(APPEND DEAL_II_COMPONENTS DOCUMENTATION)

ENDIF()

OPTION(DEAL_II_COMPONENT_EXAMPLES
  "Enable configuration and installation of the example steps. This adds a COMPONENT \"examples\" to the build system."
  ON
  )
LIST(APPEND DEAL_II_COMPONENTS EXAMPLES)

OPTION(DEAL_II_COMPONENT_PACKAGE
  "Generates additional targets for packaging deal.II"
  OFF
  )
LIST(APPEND DEAL_II_COMPONENTS PACKAGE)

OPTION(DEAL_II_COMPONENT_PYTHON_BINDINGS
  "Enable configuration and installation of the python bindings. This adds a COMPONENT \"PYTHON_BINDINGS\" to the build system."
  OFF
  )
LIST(APPEND DEAL_II_COMPONENTS PYTHON_BINDINGS)

OPTION(DEAL_II_ALLOW_AUTODETECTION
  "Allow to automatically set up features by setting all undefined DEAL_II_WITH_* variables to ON or OFF"
  ON
  )

OPTION(DEAL_II_FORCE_AUTODETECTION
  "Force feature autodetection by undefining all DEAL_II_WITH_* variables prior to configure"
  OFF
  )


########################################################################
#                                                                      #
#           Configuration options for Compilation and linking:         #
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

# This is cruel, I know. But it is better to only have a known number of
# options for CMAKE_BUILD_TYPE...
IF( NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Release" AND
    NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" AND
    NOT "${CMAKE_BUILD_TYPE}" STREQUAL "DebugRelease" )
  MESSAGE(FATAL_ERROR
    "CMAKE_BUILD_TYPE does neither match Release, Debug, nor DebugRelease!"
    )
ENDIF()

#
# Configuration behaviour:
#

OPTION(DEAL_II_ALLOW_PLATFORM_INTROSPECTION
  "Allow platform introspection for CPU command sets, SSE and AVX"
  ON
  )
MARK_AS_ADVANCED(DEAL_II_ALLOW_PLATFORM_INTROSPECTION)

OPTION(DEAL_II_WITH_SIMD
  "Allow SIMD vectorization"
  ON
  )
MARK_AS_ADVANCED(DEAL_II_WITH_SIMD)

OPTION(DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS
  "Configure sensible default CFLAGS and CXXFLAGS depending on platform, compiler and build target."
  ON
  )
MARK_AS_ADVANCED(DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS)

OPTION(DEAL_II_SETUP_COVERAGE
  "Setup debug compiler flags to provide additional test coverage information. Currently only gprof is supported."
  OFF
  )
MARK_AS_ADVANCED(DEAL_II_SETUP_COVERAGE)

OPTION(DEAL_II_UNITY_BUILD
  "Compile the library by concatenating together source files to form a few large targets instead of many small ones. This lowers total compilation wall time by about 25%."
  OFF)
MARK_AS_ADVANCED(DEAL_II_UNITY_BUILD)

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

SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH "ON" CACHE BOOL
  "Set the rpath of the library to the external link paths on installation"
  )
MARK_AS_ADVANCED(CMAKE_INSTALL_RPATH_USE_LINK_PATH)


########################################################################
#                                                                      #
#                       Compilation and linking:                       #
#                                                                      #
########################################################################

#
# Hide all unused CMake variables:
#

SET(DEAL_II_REMOVED_FLAGS
  CMAKE_CXX_FLAGS
  CMAKE_CXX_FLAGS_RELEASE
  CMAKE_CXX_FLAGS_DEBUG
  CMAKE_CXX_FLAGS_MINSIZEREL
  CMAKE_CXX_FLAGS_RELWITHDEBINFO
  CMAKE_C_FLAGS
  CMAKE_C_FLAGS_RELEASE
  CMAKE_C_FLAGS_DEBUG
  CMAKE_C_FLAGS_MINSIZEREL
  CMAKE_C_FLAGS_RELWITHDEBINFO
  CMAKE_Fortran_FLAGS
  CMAKE_Fortran_FLAGS_RELEASE
  CMAKE_Fortran_FLAGS_DEBUG
  CMAKE_Fortran_FLAGS_MINSIZEREL
  CMAKE_Fortran_FLAGS_RELWITHDEBINFO
  CUDA_NVCC_FLAGS
  CMAKE_CUDA_FLAGS
  CMAKE_CUDA_FLAGS_RELEASE
  CMAKE_CUDA_FLAGS_DEBUG
  CMAKE_CUDA_FLAGS_MINSIZEREL
  CMAKE_CUDA_FLAGS_RELWITHDEBINFO
  CMAKE_SHARED_LINKER_FLAGS
  CMAKE_SHARED_LINKER_FLAGS_DEBUG
  CMAKE_SHARED_LINKER_FLAGS_MINSIZEREL
  CMAKE_SHARED_LINKER_FLAGS_RELEASE
  CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO
  )
FOREACH(_flag ${DEAL_II_REMOVED_FLAGS})
  #
  # Promote all variables to internal cache. This prevents CMake from
  # populating these variables with default values. Further, store the
  # actual content of the variables such that users can still use
  # CMAKE_CXX_FLAGS(|_RELEASE|_DEBUG).
  #
  SET(${_flag} ${${_flag}} CACHE INTERNAL "" FORCE)
ENDFOREACH()

#
# Promote our configuration variables to cache:
#

SET(DEAL_II_USED_FLAGS
  DEAL_II_CUDA_FLAGS
  DEAL_II_CUDA_FLAGS_DEBUG
  DEAL_II_CUDA_FLAGS_RELEASE
  DEAL_II_CXX_FLAGS
  DEAL_II_CXX_FLAGS_DEBUG
  DEAL_II_CXX_FLAGS_RELEASE
  DEAL_II_LINKER_FLAGS
  DEAL_II_LINKER_FLAGS_DEBUG
  DEAL_II_LINKER_FLAGS_RELEASE
  )
FOREACH(_flag ${DEAL_II_USED_FLAGS})
  SET(${_flag} "${${_flag}}" CACHE STRING
    "The user supplied cache variable will be appended _at the end_ of the configuration step to the auto generated ${_flag} variable"
    )
  MARK_AS_ADVANCED(${_flag})
ENDFOREACH()

FOREACH(_variable
  DEAL_II_DEFINITIONS
  DEAL_II_DEFINITIONS_DEBUG
  DEAL_II_DEFINITIONS_RELEASE
  )
  SET(${_variable} ${${_variable}} CACHE STRING
    "Additional, user supplied compile definitions"
    )
  MARK_AS_ADVANCED(${_variable})
ENDFOREACH()

#
# Translate CMake specific variables to deal.II naming:
#

FOREACH(_flag
    CUDA_FLAGS CUDA_FLAGS_RELEASE CUDA_FLAGS_DEBUG
    CXX_FLAGS CXX_FLAGS_RELEASE CXX_FLAGS_DEBUG
    )
  IF(NOT "${CMAKE_${_flag}}" STREQUAL "")
    MESSAGE(STATUS
      "Prepending \${CMAKE_${_flag}} to \${DEAL_II_${_flag}}"
      )
    SET(DEAL_II_${_flag} "${CMAKE_${_flag}} ${DEAL_II_${_flag}}")
  ENDIF()
ENDFOREACH()

FOREACH(_flag LINKER_FLAGS LINKER_FLAGS_DEBUG LINKER_FLAGS_RELEASE)
  IF(NOT "${CMAKE_SHARED_${_flag}}" STREQUAL "")
    MESSAGE(STATUS
      "Prepending \${CMAKE_SHARED_${_flag}} to \${DEAL_II_${_flag}}"
      )
    SET(DEAL_II_${_flag} "${CMAKE_${_flag}} ${DEAL_II_${_flag}}")
  ENDIF()
ENDFOREACH()

IF(NOT "${CUDA_NVCC_FLAGS}}" STREQUAL "")
  MESSAGE(STATUS
    "Prepending \${CUDA_NVCC_FLAGS} to \${DEAL_II_CUDA_FLAGS}"
    )
  SET(DEAL_II_CUDA_FLAGS "${CUDA_NVCC_FLAGS} ${DEAL_II_CUDA_FLAGS}")
ENDIF()



#
# Store user supplied flags in ${_flag}_SAVED and clear configuration
# variables.
#

FOREACH(_flag ${DEAL_II_USED_FLAGS})
  #
  # The order of compiler and linker flags is important. In order to
  # provide an override mechanism we have to save the initial (cached)
  # variable at this point and clear it.
  # ${flags}_SAVED will be appended to ${flags} again in
  # setup_finalize.cmake (called at the end of the main CMakeLists.txt
  # file).
  #
  SET(${_flag}_SAVED ${${_flag}})
  SET(${_flag} "")
ENDFOREACH()

#
# Also set all unused CMAKE_* flags to an empty string for the
# configuration run so that it does not confuse the build system (to unset
# is not an option - it is cached...)
#
FOREACH(_flag ${DEAL_II_REMOVED_FLAGS})
  SET(${_flag} "")
ENDFOREACH()

#
# Finally, read in CXXFLAGS, LDFLAGS and NVCCFLAGS from environment and
# prepend them to the saved variables:
#
# Also strip leading and trailing whitespace from linker flags to make
# old cmake versions happy
#

SET(DEAL_II_CXX_FLAGS_SAVED "$ENV{CXXFLAGS} ${DEAL_II_CXX_FLAGS_SAVED}")
STRING(STRIP "${DEAL_II_CXX_FLAGS_SAVED}" DEAL_II_CXX_FLAGS_SAVED)
SET(DEAL_II_LINKER_FLAGS_SAVED "$ENV{LDFLAGS} ${DEAL_II_LINKER_FLAGS_SAVED}")
STRING(STRIP "${DEAL_II_LINKER_FLAGS_SAVED}" DEAL_II_LINKER_FLAGS_SAVED)
SET(DEAL_II_CUDA_FLAGS_SAVED "$ENV{NVCCFLAGS} ${DEAL_II_CUDA_FLAGS_SAVED}")
STRING(STRIP "${DEAL_II_CUDA_FLAGS_SAVED}" DEAL_II_CUDA_FLAGS_SAVED)
UNSET(ENV{CXXFLAGS})
UNSET(ENV{LDFLAGS})
UNSET(ENV{NVCCFLAGS})


########################################################################
#                                                                      #
#                Components and miscellaneous setup:                   #
#                                                                      #
########################################################################

OPTION(DEAL_II_WITH_64BIT_INDICES
  "If set to ON, then use 64-bit data types to represent global degree of freedom indices. The default is to OFF. You only want to set this to ON if you will solve problems with more than 2^31 (approximately 2 billion) unknowns. If set to ON, you also need to ensure that both Trilinos and/or PETSc support 64-bit indices."
  OFF
  )
LIST(APPEND DEAL_II_FEATURES 64BIT_INDICES)

OPTION(DEAL_II_WITH_COMPLEX_VALUES
  "If set to OFF, the classes that take a number type are not explicitly instantiated for std::complex<float> and std::complex<double>. This effectively disables the support for computing with complex values. If PETSc is built with complex scalar type, this option must be ON."
  ON
  )
LIST(APPEND DEAL_II_FEATURES COMPLEX_VALUES)

OPTION(DEAL_II_DOXYGEN_USE_MATHJAX
  "If set to ON, doxygen documentation is generated using mathjax"
  OFF
  )
MARK_AS_ADVANCED(DEAL_II_DOXYGEN_USE_MATHJAX)

OPTION(DEAL_II_COMPILE_EXAMPLES
  "If set to ON, all configurable example executables will be built and installed as well. If set to OFF, the examples component only installs the source code of example steps."
  ON
  )
MARK_AS_ADVANCED(DEAL_II_COMPILE_EXAMPLES)

OPTION(DEAL_II_DOXYGEN_USE_ONLINE_MATHJAX
  "If set to ON, doxygen documentation is generated using online (from CDN) mathjax copy"
  ON
  )
MARK_AS_ADVANCED(DEAL_II_DOXYGEN_USE_ONLINE_MATHJAX)

SET(DEAL_II_CPACK_EXTERNAL_LIBS "opt" CACHE STRING
    "A relative path to tree of external libraries that will be installed in bundle package. The path is relative to the /Applications/${DEAL_II_CPACK_BUNDLE_NAME}.app/Contents/Resources directory. It defaults to opt, but you may want to use a different value, for example if you want to distribute a brew based package."
  )
MARK_AS_ADVANCED(DEAL_II_CPACK_EXTERNAL_LIBS)

SET(DEAL_II_CPACK_BUNDLE_NAME "${DEAL_II_PACKAGE_NAME}" CACHE STRING
    "Name of the application bundle to generate."
  )
MARK_AS_ADVANCED(DEAL_II_CPACK_BUNDLE_NAME)


########################################################################
#                                                                      #
#                               Finalize:                              #
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
without the need to install it, if this is what you tried to do.)
"
    )
ENDIF()

#
# Miscellaneous renaming:
#

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
  IF(_var MATCHES "^(DOCUMENTATION|EXAMPLES|PACKAGE|PYTHON_BINDINGS)")
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
