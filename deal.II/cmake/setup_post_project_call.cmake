## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2013 by the deal.II authors
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


########################################################################
#                                                                      #
#         Setup that has to happen after the call to PROJECT():        #
#                                                                      #
########################################################################

#
# Library search order:
#
IF(DEAL_II_PREFER_STATIC_LIBS)
  #
  # Invert the search order for libraries when DEAL_II_PREFER_STATIC_LIBS
  # is set. This will prefer static archives instead of shared libraries:
  #
  # TODO: Does this work on a Windows or CYGWIN target?
  LIST(REVERSE CMAKE_FIND_LIBRARY_SUFFIXES)
ENDIF()


#
# Cross compilation stuff:
#
IF(CMAKE_CROSSCOMPILING)
  #
  # Disable platform introspection when cross compiling
  #
  SET(DEAL_II_ALLOW_PLATFORM_INTROSPECTION OFF CACHE BOOL "" FORCE)

  #
  # Import native expand_instantiations for use in cross compilation:
  #
  SET(DEAL_II_NATIVE "DEAL_II_NATIVE-NOTFOUND" CACHE FILEPATH
    "A pointer to a native deal.Ii build directory"
    )
  IF(DEAL_II_NATIVE MATCHES "-NOTFOUND")
    MESSAGE(FATAL_ERROR
      "Please set the CMake variable DEAL_II_NATIVE to a valid path that points to a native deal.II build directory"
      )
  ENDIF()
  INCLUDE(${DEAL_II_NATIVE}/importExecutables.cmake)
ENDIF()
