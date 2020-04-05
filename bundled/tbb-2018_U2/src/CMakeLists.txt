## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2018 by the deal.II authors
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
# Remove -Wall and -pedantic from DEAL_II_CXX_FLAGS (in directory scope) to
# avoid some annoying warnings...
#
STRIP_FLAG(DEAL_II_CXX_FLAGS "-Wall")
STRIP_FLAG(DEAL_II_CXX_FLAGS "-pedantic")

#
# As discussed in
#
#     https://software.intel.com/en-us/forums/intel-threading-building-blocks/topic/641654
#
# TBB, in a few places, will use memset() followed by placement new with the
# intent of creating objects with (at a binary level) value zero. GCC version
# 6.0 and later optimizes away the initial memset, which is not compatible with
# this approach. Hence, if supported, disable this particular optimization. The
# original TBB makefile specifies this flag for recent versions of GCC.
#
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-flifetime-dse=1")

#
# Disable other warnings:
#
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-Wimplicit-fallthrough=0")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-Wstringop-overflow=0")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-Wno-deprecated-copy")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-Wtautological-overlap-compare=0")



SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)
INCLUDE_DIRECTORIES(
  ${THREADS_BUNDLED_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/rml/include
  )

IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/version_string.ver")
  FILE(WRITE "${CMAKE_CURRENT_BINARY_DIR}/version_string.ver"
    "#define __TBB_VERSION_STRINGS(N) \"Custom\"\n"
    )
ENDIF()

SET(src_tbb
  rml/client/rml_tbb.cpp
  tbb/arena.cpp
  tbb/cache_aligned_allocator.cpp
  tbb/concurrent_hash_map.cpp
  tbb/concurrent_monitor.cpp
  tbb/concurrent_queue.cpp
  tbb/concurrent_vector.cpp
  tbb/condition_variable.cpp
  tbb/critical_section.cpp
  tbb/dynamic_link.cpp
  tbb/governor.cpp
  tbb/itt_notify.cpp
  tbb/market.cpp
  tbb/mutex.cpp
  tbb/observer_proxy.cpp
  tbb/pipeline.cpp
  tbb/private_server.cpp
  tbb/queuing_mutex.cpp
  tbb/queuing_rw_mutex.cpp
  tbb/reader_writer_lock.cpp
  tbb/recursive_mutex.cpp
  tbb/scheduler.cpp
  tbb/semaphore.cpp
  tbb/spin_mutex.cpp
  tbb/spin_rw_mutex.cpp
  tbb/task.cpp
  tbb/task_group_context.cpp
  tbb/tbb_main.cpp
  tbb/tbb_misc.cpp
  tbb/tbb_misc_ex.cpp
  tbb/tbb_statistics.cpp
  tbb/tbb_thread.cpp
  tbb/x86_rtm_rw_mutex.cpp
  )
DEAL_II_ADD_LIBRARY(obj_tbb OBJECT ${src_tbb})

#
# Add necessary definitions:
#

DEAL_II_ADD_DEFINITIONS(obj_tbb "__TBB_BUILD=1")
DEAL_II_ADD_DEFINITIONS(obj_tbb "USE_PTHREAD")
IF(NOT CMAKE_SYSTEM_NAME MATCHES "Darwin")
  DEAL_II_ADD_DEFINITIONS(obj_tbb "DO_ITT_NOTIFY")
ENDIF()
