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

SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)
INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}/rml/include
  )

IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/version_string.ver")
  FILE(WRITE "${CMAKE_CURRENT_BINARY_DIR}/version_string.ver"
    "#define __TBB_VERSION_STRINGS(N) \"Empty\"\n"
    "#define TBB_VERSION_STRINGS(N) \"Empty\"\n"
    )
ENDIF()

#
# Some architectures need additional assembly files:
#

IF(CMAKE_SYSTEM_NAME MATCHES "Windows")
  IF(CMAKE_SYSTEM_PROCESSOR MATCHES "x86")
    SET(_assembly_files
      tbb/intel64-masm/atomic_support.asm
      tbb/intel64-masm/intel64_misc.asm
      )
  ELSEIF(CMAKE_SYSTEM_PROCESSOR MATCHES "x86")
    SET(_assembly_files
      tbb/ia32-masm/atomic_support.asm
      tbb/ia32-masm/intel64_misc.asm
      )
  ENDIF()
ENDIF()

IF(CMAKE_SYSTEM_PROCESSOR MATCHES "ia64")
  SET(_assembly_files
    tbb/ia64-gas/atomic_support.s
    tbb/ia64-gas/lock_byte.s
    tbb/ia64-gas/log2.s
    tbb/ia64-gas/pause.s
    )
  IF(CMAKE_SYSTEM_NAME MATCHES "Linux")
    LIST(APPEND _assembly_files
      tbb/ia64-gas/ia64_misc.s
      )
  ENDIF()
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
  ${_assembly_files}
  )
DEAL_II_ADD_LIBRARY(obj_tbb OBJECT ${src_tbb})

#
# Add necessary definitions:
#

DEAL_II_ADD_DEFINITIONS(obj_tbb "__TBB_BUILD")

IF(CMAKE_SYSTEM_NAME MATCHES "Windows")
  DEAL_II_ADD_DEFINITIONS(obj_tbb "USE_WINTHREAD")

ELSE()

  DEAL_II_ADD_DEFINITIONS(obj_tbb "USE_PTHREAD")

  IF(NOT CMAKE_SYSTEM_NAME MATCHES "Darwin")
    DEAL_II_ADD_DEFINITIONS(obj_tbb "DO_ITT_NOTIFY")
  ENDIF()
ENDIF()

