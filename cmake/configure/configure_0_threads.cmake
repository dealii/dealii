## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2020 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------


#
# Set up threading:
#

# Clear the test flags because FindThreads.cmake will use a C compiler:
CLEAR_CMAKE_REQUIRED()

SWITCH_LIBRARY_PREFERENCE()
FIND_PACKAGE(Threads)
SWITCH_LIBRARY_PREFERENCE()

RESET_CMAKE_REQUIRED()

ADD_FLAGS(DEAL_II_LINKER_FLAGS "${CMAKE_THREAD_LIBS_INIT}")

IF(NOT Threads_FOUND)
  MESSAGE(FATAL_ERROR
    "\nFatal configuration error: CMake was unable to detect any threading "
    "support offered by the current compiler. Configuration cannot continue.\n\n"
    )
ENDIF()
