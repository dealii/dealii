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
# This file sets up:
#
#   DEAL_II_HAVE_GETHOSTNAME
#   DEAL_II_HAVE_GETPID
#   DEAL_II_HAVE_JN
#   DEAL_II_HAVE_SYS_RESOURCE_H
#   DEAL_II_HAVE_SYS_TIME_H
#   DEAL_II_HAVE_SYS_TIMES_H
#   DEAL_II_HAVE_SYS_TYPES_H
#   DEAL_II_HAVE_TIMES
#   DEAL_II_HAVE_UNISTD_H
#   DEAL_II_MSVC
#


########################################################################
#                                                                      #
#                    POSIX and Linux specific tests:                   #
#                                                                      #
########################################################################

#
# Check for various posix (and linux) specific header files and symbols
#
CHECK_INCLUDE_FILE_CXX("sys/resource.h" DEAL_II_HAVE_SYS_RESOURCE_H)

CHECK_INCLUDE_FILE_CXX("sys/time.h" DEAL_II_HAVE_SYS_TIME_H)

CHECK_INCLUDE_FILE_CXX("sys/times.h" DEAL_II_HAVE_SYS_TIMES_H)
CHECK_CXX_SYMBOL_EXISTS("times" "sys/times.h" DEAL_II_HAVE_TIMES)

CHECK_INCLUDE_FILE_CXX("sys/types.h" DEAL_II_HAVE_SYS_TYPES_H)

CHECK_INCLUDE_FILE_CXX("unistd.h" DEAL_II_HAVE_UNISTD_H)
CHECK_CXX_SYMBOL_EXISTS("gethostname" "unistd.h" DEAL_II_HAVE_GETHOSTNAME)
CHECK_CXX_SYMBOL_EXISTS("getpid" "unistd.h" DEAL_II_HAVE_GETPID)

#
# Do we have the Bessel function jn?
#
FIND_SYSTEM_LIBRARY(m_LIBRARY NAMES m)
MARK_AS_ADVANCED(m_LIBRARY)

IF(NOT m_LIBRARY MATCHES "-NOTFOUND")
  LIST(APPEND CMAKE_REQUIRED_LIBRARIES ${m_LIBRARY})
  CHECK_CXX_SYMBOL_EXISTS("jn" "math.h" DEAL_II_HAVE_JN)
  RESET_CMAKE_REQUIRED()
  IF(DEAL_II_HAVE_JN)
    LIST(APPEND DEAL_II_LIBRARIES ${m_LIBRARY})
  ENDIF()
ENDIF()


########################################################################
#                                                                      #
#                        Mac OSX specific setup:                       #
#                                                                      #
########################################################################

IF(CMAKE_SYSTEM_NAME MATCHES "Darwin")
  #
  # Use -Wno-long-double on Apple Darwin to avoid some unnecessary
  # warnings. However, newer gccs on that platform do not have
  # this flag any more, so check whether we can indeed do this
  #
  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-Wno-long-double")

  #
  # On Mac OS X, -rdynamic is accepted by the compiler (i.e.
  # it doesn't produce an error) but we always get a warning
  # that it isn't supported.
  #
  # TODO: MM: Check whether this is still necessary...
  #
  STRIP_FLAG(DEAL_II_LINKER_FLAGS "-rdynamic")

  #
  # At least on Clang 5.0.0 the template depth is set to 128, which is too low
  # to compile parts of the library. Fix this by setting a large value.
  #
  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-ftemplate-depth=1024")
ENDIF()



########################################################################
#                                                                      #
#                   Windows and CYGWIN specific setup:                 #
#                                                                      #
########################################################################


IF( CMAKE_SYSTEM_NAME MATCHES "CYGWIN"
    OR CMAKE_SYSTEM_NAME MATCHES "Windows" )
  IF( CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND
      CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8" )

    #
    # Print a warning if we have a gcc port that is older than gcc-4.8:
    #

    MESSAGE(WARNING "\n"
      "GCC ports (Cygwin, MinGW-w64, or MinGW) older than version 4.8 are unsupported on Windows\n\n"
      )

  ELSE()

    #
    # Workaround for a compiler bug on Windows platforms with modern gcc:
    #
    # GCC seems to have a hard time with
    #
    #   next_free_pair_object and next_free_single_object
    #
    # defined in tria_objects.h. It might explode in one ot the following ways:
    #  a) Internal compiler error
    #  b) emition of a truncated external symbol
    #  c) "template <...>::dimension could not be converted to 'int'"
    #
    # It seems to help to specify "-ggdb" also for optimized mode.
    #
    # TODO: Track down bug and fix properly (just kidding).
    #
    # Maier, 2013
    #

    ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_RELEASE "-g")
    ENABLE_IF_SUPPORTED(DEAL_II_LINKER_FLAGS_RELEASE "-g")
    REPLACE_FLAG(DEAL_II_CXX_FLAGS_DEBUG "-Og" "-O1")
    REPLACE_FLAG(DEAL_II_CXX_FLAGS_DEBUG "-O0" "-O1")
    REPLACE_FLAG(DEAL_II_CXX_FLAGS_DEBUG "-ggdb" "-g")
    REPLACE_FLAG(DEAL_II_LINKER_FALGS_DEBUG "-ggdb" "-g")

  ENDIF()
ENDIF()


IF(CMAKE_SYSTEM_NAME MATCHES "Windows")

  #
  # Export DEAL_II_MSVC if we are on a Windows platform:
  #
  SET(DEAL_II_MSVC TRUE)

  #
  # Shared library handling:
  #

  IF(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    # With MinGW we're lucky:
    ENABLE_IF_LINKS(DEAL_II_LINKER_FLAGS "-Wl,--export-all-symbols")
    ENABLE_IF_LINKS(DEAL_II_LINKER_FLAGS "-Wl,--enable-auto-import")
    ENABLE_IF_LINKS(DEAL_II_LINKER_FLAGS "-Wl,--allow-multiple-definition")
  ELSE()
    # Otherwise disable shared libraries:
    MESSAGE(WARNING "\n"
      "BUILD_SHARED_LIBS forced to OFF\n\n"
      )
    SET(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
  ENDIF()

ENDIF()
