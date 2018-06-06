## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2018 by the deal.II authors
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

# Notice how intelligent the version numbering of "Microsoft Visual Studio 2017
# version 15.0" is, the c++ compiler is advertised as "MSVC++ 14.1" but the
# version information is 19.10.x (this is the numbering used by CMake), see
# https://en.wikipedia.org/wiki/Microsoft_Visual_C%2B%2B#Internal_version_numbering
IF(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "19.10" )
  MESSAGE(FATAL_ERROR "\n"
    "You're using an old version of the Visual Studio C++ Compiler!\n"
    "You need at least version Visual Studio 2017.\n"
    )
ENDIF()


# enable exception handling:
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "/EHsc")


# Globally disable some legacy min and max macros that cause problems:
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "/NOMINMAX")
LIST(APPEND DEAL_II_DEFINITIONS "NOMINMAX")
LIST(APPEND DEAL_II_USER_DEFINITIONS "NOMINMAX")

# fix "fatal error C1128: number of sections exceeded object file format limit"
# happening in debug mode with visual studio 2015
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "/bigobj")

ADD_FLAGS(DEAL_II_CXX_FLAGS "/W3")

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
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "/wd4068")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "/wd4244")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "/wd4267")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "/wd4996")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "/wd4355")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "/wd4800")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "/wd4146")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "/wd4667")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "/wd4520")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "/wd4700")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "/wd4789")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "/wd4808")


# Bug in MSVC 2017: bogus warning C5037: an out-of-line definition of a member of a class template cannot have default arguments
# see https://developercommunity.visualstudio.com/content/problem/81223/incorrect-error-c5037-with-permissive.html
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "/wd5037")

#############################
#                           #
#    For Release target:    #
#                           #
#############################

IF (CMAKE_BUILD_TYPE MATCHES "Release")
  #
  # General optimization flags: (very basic for now)
  #
  ADD_FLAGS(DEAL_II_CXX_FLAGS_RELEASE "/O2")
ENDIF()


###########################
#                         #
#    For Debug target:    #
#                         #
###########################

IF (CMAKE_BUILD_TYPE MATCHES "Debug")
  LIST(APPEND DEAL_II_DEFINITIONS_DEBUG "DEBUG")
  LIST(APPEND DEAL_II_USER_DEFINITIONS_DEBUG "DEBUG")

  # generate some debug info:
  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_DEBUG "/Zi /MDd /Od")
ENDIF()

