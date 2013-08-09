## ---------------------------------------------------------------------
## $Id$
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
# This file sets up:
#
#   HAVE_GETHOSTNAME
#   HAVE_GETPID
#   HAVE_JN
#   HAVE_RAND_R
#   HAVE_SYS_RESOURCE_H
#   HAVE_SYS_TIME_H
#   HAVE_SYS_TIMES_H
#   HAVE_SYS_TYPES_H
#   HAVE_TIMES
#   HAVE_UNISTD_H
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
CHECK_INCLUDE_FILE_CXX("sys/resource.h" HAVE_SYS_RESOURCE_H)

CHECK_INCLUDE_FILE_CXX("sys/time.h" HAVE_SYS_TIME_H)

CHECK_INCLUDE_FILE_CXX("sys/times.h" HAVE_SYS_TIMES_H)
CHECK_CXX_SYMBOL_EXISTS("times" "sys/times.h" HAVE_TIMES)

CHECK_INCLUDE_FILE_CXX("sys/types.h" HAVE_SYS_TYPES_H)

CHECK_INCLUDE_FILE_CXX("unistd.h" HAVE_UNISTD_H)
CHECK_CXX_SYMBOL_EXISTS("gethostname" "unistd.h" HAVE_GETHOSTNAME)
CHECK_CXX_SYMBOL_EXISTS("getpid" "unistd.h" HAVE_GETPID)

CHECK_CXX_SYMBOL_EXISTS("rand_r" "stdlib.h" HAVE_RAND_R)

#
# Do we have the Bessel function jn?
#
# Switch the library preference back to prefer dynamic libraries if
# DEAL_II_PREFER_STATIC_LIBS=TRUE but DEAL_II_STATIC_EXECUTABLE=FALSE. In
# this case system libraries should be linked dynamically.
#
SWITCH_LIBRARY_PREFERENCE()
FIND_LIBRARY(m_LIBRARY NAMES m)
SWITCH_LIBRARY_PREFERENCE()
MARK_AS_ADVANCED(m_LIBRARY)

IF(NOT m_LIBRARY MATCHES "-NOTFOUND")
  LIST(APPEND CMAKE_REQUIRED_LIBRARIES ${m_LIBRARY})
  CHECK_CXX_SYMBOL_EXISTS("jn" "math.h" HAVE_JN)
  RESET_CMAKE_REQUIRED()
  IF(HAVE_JN)
    LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES ${m_LIBRARY})
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
  ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-Wno-long-double")

  #
  # On Mac OS X, -rdynamic is accepted by the compiler (i.e.
  # it doesn't produce an error) but we always get a warning
  # that it isn't supported.
  #
  # TODO: MM: Check whether this is still necessary...
  #
  STRIP_FLAG(DEAL_II_LINKER_FLAGS "-rdynamic")
ENDIF()



########################################################################
#                                                                      #
#                   Windows and CYGWIN specific setup:                 #
#                                                                      #
########################################################################

IF(CMAKE_SYSTEM_NAME MATCHES "Windows")
  #
  # Export DEAL_II_MSVC if we are on a Windows platform:
  #
  SET(DEAL_II_MSVC TRUE)

  #
  # Shared library handling:
  #
  IF(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    #
    # With MinGW we're lucky:
    #
    ENABLE_IF_LINKS(DEAL_II_LINKER_FLAGS "-Wl,--export-all-symbols")
    ENABLE_IF_LINKS(DEAL_II_LINKER_FLAGS "-Wl,--enable-auto-import")
    ENABLE_IF_LINKS(DEAL_II_LINKER_FLAGS "-Wl,--allow-multiple-definition")

    #
    # Workaround for a miscompilation and linkage issue with shared libraries
    # with MinGW. Replacing -O0 with -O1 seems to help..
    #
    REPLACE_FLAG(DEAL_II_CXX_FLAGS_DEBUG "-O0" "-O1")

  ELSE()

    #
    # Otherwise disable shared libraries:
    #
    MESSAGE(WARNING "\n"
      "BUILD_SHARED_LIBS forced to OFF\n\n"
      )
    SET(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
  ENDIF()

  #
  # Disable -ggdb and -g on Windows/MinGW targets for the moment until the
  # compilation issues with too big files are resolved
  #
  # - Matthias Maier, 2012
  #
  STRIP_FLAG(DEAL_II_CXX_FLAGS_DEBUG "-ggdb")
  STRIP_FLAG(DEAL_II_LINKER_FLAGS_DEBUG "-ggdb")
  STRIP_FLAG(DEAL_II_CXX_FLAGS_DEBUG "-g")
  STRIP_FLAG(DEAL_II_LINKER_FLAGS_DEBUG "-g")
ENDIF()


IF(CMAKE_SYSTEM_NAME MATCHES "CYGWIN")
  #
  # Workaround for a miscompilation and linkage issue with shared libraries
  # under Cygwin. Replacing -O0 with -O1 helps.
  #
  # - Matthias Maier, 2013
  #
  REPLACE_FLAG(DEAL_II_CXX_FLAGS_DEBUG "-O0" "-O1")
ENDIF()
