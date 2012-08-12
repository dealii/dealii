FIND_PACKAGE(Threads REQUIRED)

IF(DEAL_II_USE_CONTRIB)
  # compile and link the contrib tbb library:
  ADD_SUBDIRECTORY(contrib/tbb)

  # This sets TBB_LIBRARY and TBB_DEBUG_LIBRARY to the full path of the
  # _installed_ library location

ELSE()
  FIND_PACKAGE(TBB REQUIRED)
ENDIF()

INCLUDE_DIRECTORIES(${TBB_INCLUDE_DIR})

IF(NOT DEAL_II_USE_CONTRIB)
  LIST(APPEND deal_ii_include_paths
    ${TBB_INCLUDE_DIR}
    )
ENDIF()

LIST(APPEND deal_ii_external_libraries
  ${TBB_LIBRARY}
  )
LIST(APPEND deal_ii_external_debug_libraries
  ${TBB_DEBUG_LIBRARY}
  ) #TODO

SET(DEAL_II_USE_MT TRUE)
SET(DEAL_II_USE_MT_POSIX TRUE)
