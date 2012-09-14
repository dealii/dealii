#
# Configuration for the tbb library:
#

INCLUDE(CheckCXXSourceCompiles)

#
# Set up genereal threading. The macro will be included in
# CONFIGURE_FEATURE_TBB_EXTERNAL/CONTRIB:
#
MACRO(SETUP_THREADING var)
  FIND_PACKAGE(Threads)

  IF(Threads_FOUND)
    SET(${var} TRUE)

    #
    # We support threading. Go on and configure the rest:
    #

    SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${CMAKE_THREAD_LIBS_INIT}")

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

  ENDIF(Threads_FOUND)
ENDMACRO()



#
# Set up the tbb library:
#

MACRO(FIND_FEATURE_TBB_EXTERNAL var)

  FIND_PACKAGE(TBB)

  # In case we don't have a debug library:
  IF(NOT TBB_DEBUG_FOUND)
    SET(TBB_DEBUG_LIBRARY ${TBB_LIBRARY})
  ENDIF()

  IF(TBB_FOUND)
    #SET(${var} TRUE)
  ENDIF()

ENDMACRO()


MACRO(CONFIGURE_FEATURE_TBB_EXTERNAL var)

  INCLUDE_DIRECTORIES(${TBB_INCLUDE_DIR})

  IF (CMAKE_BUILD_TYPE MATCHES "Debug")
    LIST(APPEND deal_ii_external_libraries ${TBB_DEBUG_LIBRARY})
  ELSE()
    LIST(APPEND deal_ii_external_libraries ${TBB_LIBRARY})
  ENDIF()

  # Setup threading and if successfull return TRUE:
  SETUP_THREADING(${var})
ENDMACRO()


SET(HAVE_CONTRIB_FEATURE_TBB TRUE)


MACRO(CONFIGURE_FEATURE_TBB_CONTRIB var)
  #
  # Add tbb directly to the object files of deal.II
  #

  # Setup threading (before configurating our build...)
  # and if successfull return TRUE:
  SETUP_THREADING(${var})

  INCLUDE_DIRECTORIES(
    ${CMAKE_SOURCE_DIR}/contrib/tbb/tbb30_104oss/include
    )

  ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/contrib/tbb/tbb30_104oss/src)

  LIST(APPEND deal_ii_additional_object_files
    $<TARGET_OBJECTS:obj_tbb>
    )
ENDMACRO()


MACRO(CONFIGURE_FEATURE_TBB_ERROR_MESSAGE)
    MESSAGE(SEND_ERROR "
Could not find the tbb library!

Please ensure that the tbb library is installed on your computer.
If the library is not at a default location, either provide some hints
via environment variables:
TBB_LIBRARY_DIR TBB_INCLUDE_DIR
Or set the relevant variables by hand in ccmake.

Alternatively you may choose to compile the bundled contrib library of
boost by setting DEAL_II_ALLOW_CONTRIB=on or
DEAL_II_FORCE_CONTRIB_TBB=on.

")
ENDMACRO()


CONFIGURE_FEATURE(TBB)
