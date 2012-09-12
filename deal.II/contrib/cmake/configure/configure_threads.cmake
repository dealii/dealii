INCLUDE(CheckCXXSourceCompiles)


#
# Set up genereal threading:
#

FIND_PACKAGE(Threads REQUIRED)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_THREAD_LIBS_INIT}")
LIST(APPEND deal_ii_required_linker_flags ${CMAKE_THREAD_LIBS_INIT})

SET(DEAL_II_USE_MT TRUE)


#
# Set up some posix threads specific configuration toggles:
#
IF(CMAKE_USE_PTHREADS_INIT)
  SET(DEAL_II_USE_MT_POSIX TRUE)

  # Check whether posix thread barriers are available:

  LIST(APPEND CMAKE_REQUIRED_FLAGS "${CMAKE_THREAD_LIBS_INIT}")

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

  LIST(REMOVE_ITEM CMAKE_REQUIRED_FLAGS "${CMAKE_THREAD_LIBS_INIT}")

  IF(NOT DEAL_II_HAVE_MT_POSIX_BARRIERS)
    SET(DEAL_II_USE_MT_POSIX_NO_BARRIERS TRUE)
  ENDIF()
ENDIF()


#
# In some cases, -threads (or whatever else command line option)
# switches on some preprocessor flags. If this is not the case,
# then define them explicitely.
#

LIST(APPEND CMAKE_REQUIRED_FLAGS "${CMAKE_THREAD_LIBS_INIT}")
CHECK_CXX_SOURCE_COMPILES(
  "
  #if !defined (_REENTRANT) && !defined (_THREAD_SAFE)
  # error Neither _REENTRANT nor _THREAD_SAFE were defined.
    nonsense
  #endif
  int main(){ return 0; }
  "
  DEAL_II_HAVE_SANE_MT_DEFINITIONS)
LIST(REMOVE_ITEM CMAKE_REQUIRED_FLAGS "${CMAKE_THREAD_LIBS_INIT}")

IF(NOT DEAL_II_HAVE_SANE_MT_DEFINITIONS)
  ADD_DEFINITIONS("-D_REENTRANT" "-D_THREAD_SAFE")
ENDIF()


#
# Set up the tbb library:
#


IF(NOT DEAL_II_FORCE_CONTRIB_TBB)
  IF(DEAL_II_ALLOW_CONTRIB)
    FIND_PACKAGE(TBB)
  ELSE()
    FIND_PACKAGE(TBB REQUIRED)
  ENDIF()

  # In case we don't have a debug library:
  IF(NOT TBB_DEBUG_FOUND)
    SET(TBB_DEBUG_LIBRARY ${TBB_LIBRARY})
  ENDIF()
ENDIF()

IF(DEAL_II_FORCE_CONTRIB_TBB OR NOT TBB_FOUND)

  # TODO: The same as with everything else...

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

IF (CMAKE_BUILD_TYPE MATCHES "debug")
  LIST(APPEND deal_ii_external_libraries ${TBB_DEBUG_LIBRARY})
ELSE()
  LIST(APPEND deal_ii_external_libraries ${TBB_LIBRARY})
ENDIF()
