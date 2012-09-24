#
# Configuration for tbb support:
#


#
# Set up genereal threading:
# The macro will be included in CONFIGURE_FEATURE_TBB_EXTERNAL/CONTRIB.
#
MACRO(SETUP_THREADING var)
  FIND_PACKAGE(Threads)

  IF(Threads_FOUND)
    MARK_AS_ADVANCED(
      pthread_LIBRARY
      )
    SET(${var} TRUE)

    #
    # We support threading. Go on and configure the rest:
    #

    ADD_FLAGS(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_THREAD_LIBS_INIT}")
    SET(DEAL_II_USE_MT TRUE)

    #
    # Set up some posix thread specific configuration toggles:
    #
    IF(CMAKE_USE_PTHREADS_INIT)
      SET(DEAL_II_USE_MT_POSIX TRUE)

      #
      # Check whether posix thread barriers are available:
      #
      ADD_FLAGS(CMAKE_REQUIRED_FLAGS "${CMAKE_THREAD_LIBS_INIT}")
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
      STRIP_FLAG(CMAKE_REQUIRED_FLAGS "${CMAKE_THREAD_LIBS_INIT}")

      IF(NOT DEAL_II_HAVE_MT_POSIX_BARRIERS)
        SET(DEAL_II_USE_MT_POSIX_NO_BARRIERS TRUE)
      ENDIF()

    ENDIF()

    #
    # In some cases, -threads (or whatever else command line option)
    # switches on some preprocessor flags. If this is not the case,
    # then define them explicitely.
    #
    ADD_FLAGS(CMAKE_REQUIRED_FLAGS "${CMAKE_THREAD_LIBS_INIT}")
    CHECK_CXX_SOURCE_COMPILES(
      "
      #if !defined (_REENTRANT) && !defined (_THREAD_SAFE)
      # error Neither _REENTRANT nor _THREAD_SAFE were defined.
        nonsense
      #endif
      int main(){ return 0; }
      "
      DEAL_II_HAVE_SANE_MT_DEFINITIONS)
    STRIP_FLAG(CMAKE_REQUIRED_FLAGS "${CMAKE_THREAD_LIBS_INIT}")

    IF(NOT DEAL_II_HAVE_SANE_MT_DEFINITIONS)
      LIST(APPEND DEAL_II_DEFINITIONS "_REENTRANT" "_THREAD_SAFE")
      LIST(APPEND DEAL_II_USER_DEFINITIONS "_REENTRANT" "_THREAD_SAFE") # TODO: Necessary?
    ENDIF()

  ENDIF(Threads_FOUND)
ENDMACRO()



#
# Set up the tbb library:
#

MACRO(FEATURE_TBB_FIND_EXTERNAL var)
  FIND_PACKAGE(TBB)

  IF(TBB_FOUND)
    SET(${var} TRUE)
  ENDIF()
ENDMACRO()


MACRO(FEATURE_TBB_CONFIGURE_EXTERNAL var)

  INCLUDE_DIRECTORIES(${TBB_INCLUDE_DIR})

  IF (CMAKE_BUILD_TYPE MATCHES "Debug")
    IF(TBB_DEBUG_FOUND)
      LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES_DEBUG ${TBB_DEBUG_LIBRARY})
      LIST(APPEND DEAL_II_DEFINITIONS_DEBUG
        "TBB_USE_DEBUG=1" "TBB_DO_ASSERT=1"
        )
    ELSE()
      MESSAGE(STATUS
        "No debug tbb library was found. "
        "The regular tbb lib will be used for the debug target instead."
        )
      LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES_DEBUG ${TBB_LIBRARY})
    ENDIF()
  ENDIF()

  IF (CMAKE_BUILD_TYPE MATCHES "Release")
    LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES_RELEASE ${TBB_LIBRARY})
  ENDIF()

  # Setup threading and if successfull return TRUE:
  SETUP_THREADING(${var})
ENDMACRO()


SET(FEATURE_TBB_HAVE_CONTRIB TRUE)


MACRO(FEATURE_TBB_CONFIGURE_CONTRIB var)
  #
  # Add tbb directly to the object files of deal.II
  #

  #
  # Setup threading (before configuring our build...)
  # and if successfull return TRUE:
  #
  SETUP_THREADING(${var})

  #
  # We have to disable a bunch of warnings:
  #
  ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-Wno-parentheses")
  ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-Wno-long-long")

  #
  # Add some definitions to use the header files in debug mode:
  #
  IF (CMAKE_BUILD_TYPE MATCHES "Debug")
    LIST(APPEND DEAL_II_DEFINITIONS_DEBUG
      "TBB_DO_DEBUG=1" "TBB_DO_ASSERT=1"
      )
  ENDIF()

  INCLUDE_DIRECTORIES(
    ${CMAKE_SOURCE_DIR}/contrib/tbb/tbb30_104oss/include
    )

  ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/contrib/tbb/tbb30_104oss/src)

ENDMACRO()


CONFIGURE_FEATURE(TBB)

