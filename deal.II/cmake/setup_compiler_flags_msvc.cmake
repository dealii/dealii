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
## Author: Timo Heister
##
#####

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
ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "/EHsc")

#enable warnings:
ADD_FLAGS(CMAKE_CXX_FLAGS "/W3")

# Globally disable some legacy min and max macros that cause problems:
ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "/NOMINMAX")

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


#########################
#                       #
#    Set up C FLAGS:    #
#                       #
#########################

#
# For the moment we assume that CC and CXX are the same compiler and that
# we can set (almost) the same default flags for both:
#
SET(CMAKE_C_FLAGS ${CMAKE_CXX_FLAGS})
SET(DEAL_II_C_FLAGS_RELEASE ${DEAL_II_CXX_FLAGS_RELEASE})
SET(DEAL_II_C_FLAGS_DEBUG ${DEAL_II_CXX_FLAGS_DEBUG})

