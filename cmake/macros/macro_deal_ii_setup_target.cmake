## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2012 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

#
# This file implements the DEAL_II_SETUP_TARGET macro, which is
# part of the deal.II library.
#
# Usage:
#       deal_ii_setup_target(target)
#       deal_ii_setup_target(target [DEBUG|RELEASE] [MODULE])
#
# This macro appends the deal.II target to the link interface of the specified
# target, which in turn ensures that necessary include directories, linker
# flags, compile flags and compile definitions are set.
#
# DEBUG or RELEASE build type
# ---------------------------
# If no "DEBUG" or "RELEASE" keyword is specified after the target, the
# current CMAKE_BUILD_TYPE is used instead. A CMAKE_BUILD_TYPE "Debug" is
# equivalent to the DEBUG keyword, a CMAKE_BUILD_TYPE "Release" is
# equivalent to the RELEASE keyword.
#
# This macro throws a FATAL_ERROR in case no DEBUG/RELEASE keyword is set
# and the build type is different from "Debug", or "Release".
#
# If the requested build type is not available (e.g. DEBUG request but
# deal.II was compiled with release mode only), the macro also throws a
# FATAL_ERROR.
#
# MODULE build type
# -----------------
# If the MODULE argument is provided, then the macro links the target against
# the version of the library that is built into a C++20 module. In that case,
# the target should use
#   import dealii;
# instead of
#   #include <deal.II/...>
# to learn about deal.II's functions and classes.
#
# The macro will throw a FATAL_ERROR if MODULE is specified among the arguments,
# but the library was not built with modules.
#

macro(deal_ii_setup_target _target)

  if(NOT DEAL_II_PROJECT_CONFIG_INCLUDED)
    message(FATAL_ERROR
      "\nDEAL_II_SETUP_TARGET can only be called in external projects after "
      "the inclusion of deal.IIConfig.cmake. It is not intended for "
      "internal use.\n\n"
      )
  endif()

  if(NOT DEAL_II_TARGET_CONFIG_INCLUDED)
    include(${DEAL_II_TARGET_CONFIG})
    set(DEAL_II_TARGET_CONFIG_INCLUDED TRUE)
  endif()

  #
  # Parse the optional arguments passed to this macro. These
  # may only contain the build type (DEBUG or RELEASE), and the
  # MODULE keyword.
  #
  set(_build)
  set(_module)
  foreach (_arg ${ARGN})
    if("${_arg}" MATCHES "^(DEBUG|RELEASE)$")
      if("${_build}" STREQUAL "")
        set(_build ${_arg})
      else()
        message(FATAL_ERROR
          "\nThe deal_ii_setup_target() macro can take optional arguments, "
          "but the debug or release configuration can only be specified "
          "once. The arguments you passed are \"${ARGN}\"."
          "\n\n"
          )
      endif()
    elseif("${_arg}" STREQUAL "MODULE")
      if("${_module}" STREQUAL "")
        if(DEAL_II_BUILD_CXX20_MODULE)
          set(_module "ON")
        else()
          message(FATAL_ERROR
            "\nYou provided the MODULE option to the deal_ii_setup_target() macro, "
            "but this option can only be used if the library was actually built "
            "with C++20 modules enabled."
            "\n\n"
            )
        endif()
      else()
        message(FATAL_ERROR
          "\nThe deal_ii_setup_target() macro can take optional arguments, "
          "but the MODULE option can only be specified once. The arguments "
          "you passed are <${ARGN}>."
          "\n\n"
          )
      endif()
    else()
      message(FATAL_ERROR
        "\nThe deal_ii_setup_target() macro was called with an invalid argument. "
        "Valid arguments are (none), DEBUG or RELEASE, and MODULE. "
        "The argument given is \"${_arg}\"."
        "\n\n"
        )
    endif()
  endforeach()

  #
  # Set build type if not already set above by explicitly specifying
  # it as an additional macro argument that must be either DEBUG or
  # RELEASE (parsed above); alternatively, failing a specific optional
  # argument, it is set based on CMAKE_BUILD_TYPE:
  #
  if("${_build}" STREQUAL "")
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
      set(_build "DEBUG")
    elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
      set(_build "RELEASE")
    else()
      message(FATAL_ERROR
        "\nDEAL_II_SETUP_TARGET cannot determine DEBUG, or RELEASE flavor "
        "for target. CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\" is neither "
        "equal to \"Debug\", nor \"Release\"\n"
        "Set CMAKE_BUILD_TYPE accordingly, or use an explicit annotation: "
        "  deal_ii_setup_target(<target> DEBUG|RELEASE)\n\n"
        )
    endif()
  endif()

  #
  # We can only append the debug or release interface if deal.II was built
  # with the Debug or DebugRelease build type. So test for this:
  #

  if("${_build}" STREQUAL "DEBUG" AND NOT DEAL_II_BUILD_TYPE MATCHES "Debug")
    set(_build "RELEASE")
  endif()

  target_compile_flags(${_target} PRIVATE "$<COMPILE_LANGUAGE:CXX>"
    "${DEAL_II_WARNING_FLAGS} ${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${_build}}"
    )

  get_property(_type TARGET ${_target} PROPERTY TYPE)
  if(NOT "${_type}" STREQUAL "OBJECT_LIBRARY")
    target_link_flags(${_target} PRIVATE
      "${DEAL_II_LINKER_FLAGS} ${DEAL_II_LINKER_FLAGS_${_build}}"
      )
  endif()

  #
  # Now link against the right deal.II library target
  #
  if("${_module}" STREQUAL "")
    target_link_libraries(${_target} ${DEAL_II_TARGET_${_build}})
  else()
    target_link_libraries(${_target} ${DEAL_II_TARGET_NAME}_module_${_build}})
    # Make sure CMake figures out that we need to scan the sources for
    # module. This should happen automatically via policy 0155, but
    # doesn't appear to work with CMake 3.28. In any case, it doesn't
    # hurt to be explicit.
    set_property(TARGET ${_target} PROPERTY CXX_SCAN_FOR_MODULES ON)
  endif()
endmacro()
