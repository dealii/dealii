FIND_PACKAGE(Threads REQUIRED)

IF(DEAL_II_USE_CONTRIB)
  # Compiles and links libtbb, exports the TBB_* variables as well:
  ADD_SUBDIRECTORY(contrib/tbb)
ELSE()
  FIND_PACKAGE(TBB REQUIRED)
ENDIF()

INCLUDE_DIRECTORIES(${TBB_INCLUDE_DIR})

SET(deal_ii_external_libraries
  ${deal_ii_external_libraries}
  ${TBB_LIBRARY}
  )
SET(deal_ii_external_debug_libraries
  ${deal_ii_external_debug_libraries}
  ${TBB_DEBUG_LIBRARY}) #TODO

SET(DEAL_II_USE_MT TRUE)
SET(DEAL_II_USE_MT_POSIX TRUE)
