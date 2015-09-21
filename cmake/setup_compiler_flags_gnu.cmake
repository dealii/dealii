## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2015 by the deal.II authors
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
# General setup for GCC and compilers sufficiently close to GCC
#
# Please read the fat note in setup_compiler_flags.cmake prior to
# editing this file.
#

IF( CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND
    CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.4" )
  MESSAGE(WARNING "\n"
    "You're using an old version of the GNU Compiler Collection (gcc/g++)!\n"
    "It is strongly recommended to use at least version 3.4.\n"
    )
ENDIF()


########################
#                      #
#    General setup:    #
#                      #
########################

#
# Set -pedantic if the compiler supports it.
#
IF(NOT (CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND
        CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.4"))
  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-pedantic")
ENDIF()

#
# Set the pic flag.
#
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-fpic")

#
# Check whether the -as-needed flag is available. If so set it to link
# the deal.II library with it.
#
ENABLE_IF_LINKS(DEAL_II_LINKER_FLAGS "-Wl,--as-needed")

#
# Setup various warnings:
#
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-Wall")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-Wextra")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-Wpointer-arith")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-Wwrite-strings")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-Wsynth")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-Wsign-compare")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-Wswitch")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-Woverloaded-virtual")

#
# Disable Wlong-long that will trigger a lot of warnings when compiling
# with disabled C++11 support:
#
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-Wno-long-long")

#
# Disable deprecation warnings
#
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-Wno-deprecated-declarations")

IF(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  #
  # Silence Clang warnings about unused compiler parameters (works around a
  # regression in the clang driver frontend of certain versions):
  #
  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-Qunused-arguments")

  #
  # Clang verbosely warns about not supporting all our friend declarations
  # (and consequently removing access control altogether)
  #
  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-Wno-unsupported-friend")

  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-Wno-unused-parameter")
  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-Wno-unused-variable")

  # without c++11 enabled, clang produces a ton of warnings in boost:
  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-Wno-c99-extensions")
  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-Wno-variadic-macros")
  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-Wno-c++11-extensions")

  #
  # Clang versions prior to 3.6 emit a lot of false positives wrt
  # "-Wunused-function". Also suppress warnings for Xcode older than 6.3
  # (which is equivalent to clang < 3.6).
  #
  # FIXME: I wait for the day with a clang version "4.0"... and I will
  # curse the person that thought it is a _great_ idea to come up with
  # independent version numbers for clang on Mac...
  #
  IF( CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.6" OR
      ( NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.0" AND
        CMAKE_CXX_COMPILER_VERSION VERSION_LESS "6.3") )
    ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-Wno-unused-function")
  ENDIF()
ENDIF()


IF(DEAL_II_STATIC_EXECUTABLE)
  #
  # To produce a static executable, we have to statically link libstdc++
  # and gcc's support libraries and glibc:
  #
  ENABLE_IF_SUPPORTED(DEAL_II_LINKER_FLAGS "-static")
  ENABLE_IF_SUPPORTED(DEAL_II_LINKER_FLAGS "-pthread")
ENDIF()


#############################
#                           #
#    For Release target:    #
#                           #
#############################

IF (CMAKE_BUILD_TYPE MATCHES "Release")
  #
  # General optimization flags:
  #
  ADD_FLAGS(DEAL_II_CXX_FLAGS_RELEASE "-O2")

  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_RELEASE "-funroll-loops")
  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_RELEASE "-funroll-all-loops")
  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_RELEASE "-fstrict-aliasing")
ENDIF()


###########################
#                         #
#    For Debug target:    #
#                         #
###########################

IF (CMAKE_BUILD_TYPE MATCHES "Debug")

  LIST(APPEND DEAL_II_DEFINITIONS_DEBUG "DEBUG")
  LIST(APPEND DEAL_II_USER_DEFINITIONS_DEBUG "DEBUG")

  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_DEBUG "-Og")
  #
  # If -Og is not available, fall back to -O0:
  #
  IF(NOT DEAL_II_HAVE_FLAG_Og)
    ADD_FLAGS(DEAL_II_CXX_FLAGS_DEBUG "-O0")
  ENDIF()

  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_DEBUG "-ggdb")
  ENABLE_IF_SUPPORTED(DEAL_II_LINKER_FLAGS_DEBUG "-ggdb")
  #
  # If -ggdb is not available, fall back to -g:
  #
  IF(NOT DEAL_II_HAVE_FLAG_ggdb)
    ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_DEBUG "-g")
    ENABLE_IF_SUPPORTED(DEAL_II_LINKER_FLAGS_DEBUG "-g")
  ENDIF()

  IF(DEAL_II_SETUP_COVERAGE)
    #
    # Enable test coverage
    #
    ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_DEBUG "-fno-elide-constructors")
    ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_DEBUG "-ftest-coverage -fprofile-arcs")
    ENABLE_IF_SUPPORTED(DEAL_II_LINKER_FLAGS_DEBUG "-ftest-coverage -fprofile-arcs")
  ENDIF()

ENDIF()
