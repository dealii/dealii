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

#
# Check for various posix specific header files:
#

CHECK_INCLUDE_FILE("unistd.h" HAVE_UNISTD_H)

CHECK_INCLUDE_FILE("sys/resource.h"  HAVE_SYS_RESOURCE_H)

CHECK_INCLUDE_FILE("sys/stat.h" HAVE_SYS_STAT_H)

CHECK_INCLUDE_FILE("sys/syscall.h" HAVE_SYS_SYSCALL_H)

CHECK_INCLUDE_FILE("sys/times.h" HAVE_SYS_TIMES_H)

CHECK_INCLUDE_FILE("sys/types.h" HAVE_SYS_TYPES_H)



#
# Check for various posix specific functions:
#

CHECK_FUNCTION_EXISTS(gethostname HAVE_GETHOSTNAME)

CHECK_FUNCTION_EXISTS(getpid HAVE_GETPID)

CHECK_FUNCTION_EXISTS(rand_r HAVE_RAND_R)

CHECK_FUNCTION_EXISTS(times HAVE_TIMES)

TEST_BIG_ENDIAN(DEAL_II_WORDS_BIGENDIAN)

