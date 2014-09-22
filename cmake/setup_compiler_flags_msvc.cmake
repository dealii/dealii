## ---------------------------------------------------------------------
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
# General setup for the Microsoft Visual Studio C++ Compiler (Windows)
#
# Please read the fat note in setup_compiler_flags.cmake prior to
# editing this file.
#

#TODO: this check is not working, my version is 17.0.51106.1 (2012)
#IF(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "17.0.0.0" )
##  MESSAGE(WARNING "\n"
#    "You're using an old version of the MSVC C++ Compiler!\n"
#    "It is strongly recommended to use at least version 2012.\n"
#    )
#ENDIF()


########################
#                      #
#    General setup:    #
#                      #
########################

# enable exception handling:
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "/EHsc")

#enable warnings:
ADD_FLAGS(DEAL_II_CXX_FLAGS "/W3")

# Globally disable some legacy min and max macros that cause problems:
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "/NOMINMAX")

# Disable warning about unknown pragmas
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "/wd4068")

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
  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_DEBUG "/Zi")
ENDIF()

