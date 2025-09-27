## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2013 - 2024 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

#
# This file sets up:
#
#   DEAL_II_HAVE_GETHOSTNAME
#   DEAL_II_HAVE_GETPID
#   DEAL_II_HAVE_SYS_RESOURCE_H
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

CHECK_INCLUDE_FILE_CXX("unistd.h" DEAL_II_HAVE_UNISTD_H)
CHECK_CXX_SYMBOL_EXISTS("gethostname" "unistd.h" DEAL_II_HAVE_GETHOSTNAME)
CHECK_CXX_SYMBOL_EXISTS("getpid" "unistd.h" DEAL_II_HAVE_GETPID)

########################################################################
#                                                                      #
#                        Mac OSX specific setup:                       #
#                                                                      #
########################################################################

if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
  #
  # Use -Wno-long-double on Apple Darwin to avoid some unnecessary
  # warnings. However, newer gccs on that platform do not have
  # this flag any more, so check whether we can indeed do this
  #
  enable_if_supported(DEAL_II_CXX_FLAGS "-Wno-long-double")

  #
  # On Mac OS X, -rdynamic is accepted by the compiler (i.e.
  # it doesn't produce an error) but we always get a warning
  # that it isn't supported.
  #
  # TODO: MM: Check whether this is still necessary...
  #
  strip_flag(DEAL_II_LINKER_FLAGS "-rdynamic")

  #
  # At least on Clang 5.0.0 the template depth is set to 128, which is too low
  # to compile parts of the library. Fix this by setting a large value.
  #
  enable_if_supported(DEAL_II_CXX_FLAGS "-ftemplate-depth=1024")
endif()



########################################################################
#                                                                      #
#                   Windows and CYGWIN specific setup:                 #
#                                                                      #
########################################################################

#
# Put an end to user's suffering from cygwin's defects
#
if( CMAKE_SYSTEM_NAME MATCHES "CYGWIN" OR
    CMAKE_SYSTEM_NAME MATCHES "Windows" )
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    message(FATAL_ERROR
      "\nCygwin and forks such as MinGW and MinGW-64 are unsupported due to "
      "multiple unresolved miscompilation issues.\n\n"
      )
  endif()
endif()

if(CMAKE_SYSTEM_NAME MATCHES "Windows")
  #
  # Export DEAL_II_MSVC if we are on a Windows platform:
  #
  set(DEAL_II_MSVC TRUE)

  #
  # Shared library handling:
  #
  # We disabled dynamic linking on Windows. The *.dll format only allows
  # up to 65535 objects or members, i.e., it is limited by a 16bit
  # descriptor. We need more than that for deal.II. For more information,
  # look up the linker tools error LNK1189.
  #
  # Unfortunately, this means that we are stuck with static linking.
  # As a consequence each binary will be very large.
  #
  message(WARNING "\n"
    "BUILD_SHARED_LIBS forced to OFF\n\n"
    )
  set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)

  #
  # In case we find a solution to enable dynamic linking in the future,
  # we will most probably want to use the CMake infrastructure to
  # automatically export symbols on Windows targets.
  #
  if(BUILD_SHARED_LIBS)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
  endif()
endif()
