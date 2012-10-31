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
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

#
# Check for various system features:
#

INCLUDE(TestBigEndian)
TEST_BIG_ENDIAN(DEAL_II_WORDS_BIGENDIAN)


#
# Check for various posix (and linux) specific header files:
#

CHECK_INCLUDE_FILE("sys/resource.h"  HAVE_SYS_RESOURCE_H)
CHECK_INCLUDE_FILE("sys/time.h" HAVE_SYS_TIME_H)
CHECK_INCLUDE_FILE("sys/times.h" HAVE_SYS_TIMES_H)
CHECK_INCLUDE_FILE("sys/types.h" HAVE_SYS_TYPES_H)


#
# Check for various posix specific functions. On a posix system they should
# be all defined in unistd.h. On other platforms, most notably
# Windows/MinGW unistd.h is available but not all posix functions. So test
# for each funtion as well.
#

CHECK_INCLUDE_FILE("unistd.h" HAVE_UNISTD_H)
CHECK_FUNCTION_EXISTS(gethostname HAVE_GETHOSTNAME)
CHECK_FUNCTION_EXISTS(getpid HAVE_GETPID)
CHECK_FUNCTION_EXISTS(rand_r HAVE_RAND_R)
CHECK_FUNCTION_EXISTS(times HAVE_TIMES)


#
# Export DEAL_II_MSVC if we are on a Windows platform.
#

IF(CMAKE_SYSTEM_NAME MATCHES "Windows")
  SET(DEAL_II_MSVC TRUE)
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
