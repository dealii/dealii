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

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}/rml/include
  )

IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/version_string.tmp")
  FILE(WRITE "${CMAKE_CURRENT_BINARY_DIR}/version_string.tmp"
    "#define __TBB_VERSION_STRINGS \"Empty\""
    )
ENDIF()

SET(src_tbb
  rml/client/omp_dynamic_link.cpp
  rml/client/rml_omp.cpp
  rml/client/rml_tbb.cpp
  rml/server/rml_server.cpp
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
  tbb/spin_mutex.cpp
  tbb/spin_rw_mutex.cpp
  tbb/task.cpp
  tbb/task_group_context.cpp
  tbb/tbb_main.cpp
  tbb/tbb_misc.cpp
  tbb/tbb_statistics.cpp
  tbb/tbb_thread.cpp
  )

#
# Add necessary definitions:
#
SET(obj_tbb_DEFINITIONS)
IF(NOT CMAKE_SYSTEM_NAME MATCHES "Darwin")
  LIST(APPEND obj_tbb_DEFINITIONS "DO_ITT_NOTIFY")
ENDIF()

IF(DEAL_II_USE_MT_POSIX)
  LIST(APPEND obj_tbb_DEFINITIONS "USE_PTHREAD")
ENDIF()

IF(CMAKE_SYSTEM_NAME MATCHES "CYGWIN")
  MESSAGE(FATAL_ERROR "No support for internal tbb and CYGWIN atm.")
  #LIST(APPEND obj_tbb_DEFINITIONS
  #  "USE_WINTHREAD"
  #  "_WIN32"
  #  "__MSVCRT_VERSION__=0x0700"
  #  )
  #ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-msse")
  #ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-mthreads")
ENDIF()

DEAL_II_ADD_LIBRARY(obj_tbb OBJECT ${src_tbb})

DEAL_II_ADD_DEFINITIONS(obj_tbb ${obj_tbb_DEFINITIONS})

