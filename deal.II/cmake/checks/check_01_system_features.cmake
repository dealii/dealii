#####
##
## Copyright (C) 2012, 2013 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

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


###########################################################################
#                                                                         #
#                     POSIX and Linux specific tests:                     #
#                                                                         #
###########################################################################

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
FIND_LIBRARY(m_lib NAMES m)
MARK_AS_ADVANCED(m_lib)

IF(NOT m_lib MATCHES "-NOTFOUND")
  SET(CMAKE_REQUIRED_LIBRARIES ${m_lib})
  CHECK_CXX_SYMBOL_EXISTS("jn" "math.h" HAVE_JN)
  SET(CMAKE_REQUIRED_LIBRARIES)
  IF(HAVE_JN)
    LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES ${m_lib})
  ENDIF()
ENDIF()


###########################################################################
#                                                                         #
#                    Windows and CYGWIN specific setup:                   #
#                                                                         #
###########################################################################

IF(CMAKE_SYSTEM_NAME MATCHES "Windows")
  #
  # Export DEAL_II_MSVC if we are on a Windows platform:
  #
  SET(DEAL_II_MSVC TRUE)

  #
  # Disable -ggdb and -g on Windows/MinGW targets for the moment until the
  # compilation issues with too big files are resolved
  #
  # - Matthias Maier, 2012
  #
  STRIP_FLAG(DEAL_II_CXX_FLAGS_DEBUG "-ggdb")
  STRIP_FLAG(DEAL_II_SHARED_LINKER_FLAGS_DEBUG "-ggdb")
  STRIP_FLAG(DEAL_II_CXX_FLAGS_DEBUG "-g")
  STRIP_FLAG(DEAL_II_SHARED_LINKER_FLAGS_DEBUG "-g")
ENDIF()


#
# Disable shared libraries on CYGWIN and Windows targets for the moment.
# Our support for shared libraries on Windows is a bit buggy atm..
#
# - Matthias Maier, 2012
#
IF( CMAKE_SYSTEM_NAME MATCHES "CYGWIN" OR
    CMAKE_SYSTEM_NAME MATCHES "Windows" )
  MESSAGE(WARNING "\n"
    "BUILD_SHARED_LIBS forced to OFF\n\n"
    )
  SET(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
ENDIF()

