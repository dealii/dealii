## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2022 by the deal.II authors
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
# Set up threading:
#

# Clear the test flags because FindThreads.cmake will use a C compiler:
clear_cmake_required()

#
# Ensure that we use "-pthread" instead of "-lpthread". We require
# "-pthread" for certain versions of gcc (and standard library thread
# support). Otherwise the order of libraries might be wrong on the final
# link interface and we end up with linker errors such as the following:
#   /usr/bin/ld: [...]/step-69.cc.o: undefined reference to symbol 'pthread_create@@GLIBC_2.2.5'
#   /usr/bin/ld: /usr/lib/x86_64-linux-gnu/libpthread.so: error adding symbols: DSO missing from command line
#
set(THREADS_PREFER_PTHREAD_FLAG true)

find_package(Threads)
reset_cmake_required()

add_flags(DEAL_II_LINKER_FLAGS "${CMAKE_THREAD_LIBS_INIT}")

#
# Make sure that we compile with "-pthread" as well. The "-pthread"
# compiler flag might add certain preprocessor definitions when compiling.
#
if("${CMAKE_THREAD_LIBS_INIT}" MATCHES "-pthread")
  add_flags(DEAL_II_CXX_FLAGS "${CMAKE_THREAD_LIBS_INIT}")
endif()

if(NOT Threads_FOUND)
  message(FATAL_ERROR
    "\nFatal configuration error: CMake was unable to detect any threading "
    "support offered by the current compiler. Configuration cannot continue.\n\n"
    )
endif()
