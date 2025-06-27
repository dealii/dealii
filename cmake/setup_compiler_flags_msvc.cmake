## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2013 - 2024 by the deal.II authors
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
# General setup for the Microsoft Visual Studio C++ Compiler (Windows)
#
# Please read the fat note in setup_compiler_flags.cmake prior to
# editing this file.
#


########################
#                      #
#    General setup:    #
#                      #
########################

# Notice how intelligent the version numbering of "Microsoft Visual Studio 2019
# version 16.0" is, the c++ compiler is advertised as "MSVC++ 14.20" but the
# version information is 19.20.x (this is the numbering used by CMake), see
# https://en.wikipedia.org/wiki/Microsoft_Visual_C%2B%2B#Internal_version_numbering
if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "19.20" )
  message(FATAL_ERROR "\n"
    "You're using an old version of the Visual Studio C++ Compiler!\n"
    "You need at least version Visual Studio 2019.\n"
    )
endif()


# enable exception handling:
enable_if_supported(DEAL_II_CXX_FLAGS "/EHsc")


# Globally disable some legacy min and max macros that cause problems:
enable_if_supported(DEAL_II_CXX_FLAGS "/NOMINMAX")
list(APPEND DEAL_II_DEFINITIONS "NOMINMAX")

# fix "fatal error C1128: number of sections exceeded object file format limit"
# happening in debug mode with visual studio 2015
enable_if_supported(DEAL_II_CXX_FLAGS "/bigobj")

add_flags(DEAL_II_WARNING_FLAGS "/W3")

#
# Selectively disable a bunch of warnings:
#
# 4068 - unknown pragmas
# 4244 - implied downcasting from double to float
# 4267 - implied downcasting from size_t to unsigned int
# 4996 - unsafe functions, such as strcat and sprintf
# 4355 - 'this' : used in base member initializer list
# 4661 - no suitable definition provided for explicit template instantiation request
# 4800 - forcing value to bool 'true' or 'false' (performance warning)
# 4146 - unary minus operator applied to unsigned type, result still unsigned
# 4667 - no function template defined that matches forced instantiation
# 4520 - multiple default constructors specified
# 4700 - uninitialized local variable
# 4789 - destination of memory copy is too small
# 4808 - case 'value' is not a valid value for switch condition of type 'bool
enable_if_supported(DEAL_II_WARNING_FLAGS "/wd4068")
enable_if_supported(DEAL_II_WARNING_FLAGS "/wd4244")
enable_if_supported(DEAL_II_WARNING_FLAGS "/wd4267")
enable_if_supported(DEAL_II_WARNING_FLAGS "/wd4996")
enable_if_supported(DEAL_II_WARNING_FLAGS "/wd4355")
enable_if_supported(DEAL_II_WARNING_FLAGS "/wd4800")
enable_if_supported(DEAL_II_WARNING_FLAGS "/wd4146")
enable_if_supported(DEAL_II_WARNING_FLAGS "/wd4667")
enable_if_supported(DEAL_II_WARNING_FLAGS "/wd4520")
enable_if_supported(DEAL_II_WARNING_FLAGS "/wd4700")
enable_if_supported(DEAL_II_WARNING_FLAGS "/wd4789")
enable_if_supported(DEAL_II_WARNING_FLAGS "/wd4808")


# Bug in MSVC 2017: bogus warning C5037: an out-of-line definition of a member of a class template cannot have default arguments
# see https://developercommunity.visualstudio.com/content/problem/81223/incorrect-error-c5037-with-permissive.html
enable_if_supported(DEAL_II_WARNING_FLAGS "/wd5037")

#############################
#                           #
#    For Release target:    #
#                           #
#############################

if (CMAKE_BUILD_TYPE MATCHES "Release")
  #
  # General optimization flags: (very basic for now)
  #
  add_flags(DEAL_II_CXX_FLAGS_RELEASE "/O2")

  #
  # Disable assert() in deal.II and user projects in release mode
  #
  list(APPEND DEAL_II_DEFINITIONS_RELEASE "NDEBUG")

endif()


###########################
#                         #
#    For Debug target:    #
#                         #
###########################

if (CMAKE_BUILD_TYPE MATCHES "Debug")
  list(APPEND DEAL_II_DEFINITIONS_DEBUG "DEBUG")

  # generate some debug info:
  enable_if_supported(DEAL_II_CXX_FLAGS_DEBUG "/Zi /MDd /Od")
endif()
