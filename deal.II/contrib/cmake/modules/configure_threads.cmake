FIND_PACKAGE(Threads REQUIRED)

# TODO: Necessary link commands for threads? Or are they set automatically
# by the package?

IF(DEAL_II_ALLOW_CONTRIB)
  FIND_PACKAGE(TBB)
ELSE()
  FIND_PACKAGE(TBB REQUIRED)
ENDIF()

IF(DEAL_II_FORCE_CONTRIB_TBB OR NOT TBB_FOUND)
  SET(libtbb_directory "tbb30_104oss")

  # compile and link the contrib tbb library:
  ADD_SUBDIRECTORY(contrib/tbb)


  # set TBB_LIBRARY and TBB_DEBUG_LIBRARY to the full path of the
  # _installed_ library location:

  SET(TBB_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/contrib/tbb/${libtbb_directory}/include)
  SET(TBB_LIBRARY ${CMAKE_INSTALL_PREFIX}/lib/libtbb.so)
  SET(TBB_DEBUG_LIBRARY ${CMAKE_INSTALL_PREFIX}/lib/libtbb_debug.so)
ENDIF()

INCLUDE_DIRECTORIES(${TBB_INCLUDE_DIR})

LIST(APPEND deal_ii_external_libraries
  ${TBB_LIBRARY}
  )
LIST(APPEND deal_ii_external_debug_libraries
  ${TBB_DEBUG_LIBRARY}
  )

# TODO: Renaming!
SET(DEAL_II_USE_MT TRUE)
SET(DEAL_II_USE_MT_POSIX TRUE)
