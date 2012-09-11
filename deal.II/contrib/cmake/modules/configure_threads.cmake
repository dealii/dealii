INCLUDE(CheckCXXSourceCompiles)

FIND_PACKAGE(Threads REQUIRED)

LIST(APPEND deal_ii_required_flags ${CMAKE_THREAD_LIBS_INIT})

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


SET(DEAL_II_USE_MT TRUE)

IF(CMAKE_USE_PTHREADS_INIT)
  SET(DEAL_II_USE_MT_POSIX TRUE)

  # Check whether posix thread barriers are available:

  SET(CMAKE_REQUIRED_FLAGS "${CMAKE_THREAD_LIBS_INIT}")
  CHECK_CXX_SOURCE_COMPILES(
  "
  #include <pthread.h>
  int main()
  {
    pthread_barrier_t pb;
    pthread_barrier_init (&pb, 0, 1);
    pthread_barrier_wait (&pb);
    pthread_barrier_destroy (&pb);
    return 0;
  }
  "
  DEAL_II_HAVE_MT_POSIX_BARRIERS)

  SET(CMAKE_REQUIRED_FLAGS "")

  IF(NOT DEAL_II_HAVE_MT_POSIX_BARRIERS)
    SET(DEAL_II_USE_MT_POSIX_NO_BARRIERS TRUE)
  ENDIF()

ENDIF()

