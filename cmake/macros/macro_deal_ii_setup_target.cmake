## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2025 by the deal.II authors
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
# This file implements the DEAL_II_SETUP_TARGET macro, which is
# part of the deal.II library.
#
# Usage:
#       deal_ii_setup_target(target)
#       deal_ii_setup_target(target DEBUG|RELEASE)
#
# This appends the deal.II target to the link interface of the specified
# target, which in turn ensures that necessary include directories, linker
# flags, compile flags and compile definitions are set.
#
# If no "DEBUG" or "RELEASE" keyword is specified after the target, the
# current CMAKE_BUILD_TYPE is used instead. A CMAKE_BUILD_TYPE "Debug" is
# equivalent to the DEBUG keyword, a CMAKE_BUILD_TYPE "Release" is
# equivalent to the RELEASE keyword.
#
# This macro throws a FATAL_ERROR in case no DEBUG/RELEASE keyword is set
# and the build type is different from "Debug", or "Release".
#
# If the requested build type is not available (e.g. DEBUG request but
# deal.II was compiled with release mode only), the macro throws a
# FATAL_ERROR.
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
  # may only contain the build type (DEBUG or RELEASE).
  #
  set(_build)
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
    else()
      message(FATAL_ERROR
        "\nThe deal_ii_setup_target() macro was called with an invalid argument. "
        "Valid arguments are (none), DEBUG, or RELEASE. "
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

  target_link_libraries(${_target} ${DEAL_II_TARGET_${_build}})
endmacro()
