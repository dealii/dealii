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
    SET(${var} TRUE)

    #
    # We support threading. Go on and configure the rest:
    #
    ADD_FLAGS(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_THREAD_LIBS_INIT}")
    SET(DEAL_II_USE_MT TRUE)

    #
    # Set up some posix threads specific configuration toggles:
    #
    IF(CMAKE_USE_PTHREADS_INIT)
      SET(DEAL_II_USE_MT_POSIX TRUE)

      #
      # Check whether posix thread barriers are available:
      #
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
      LIST(APPEND deal_ii_external_libraries ${TBB_DEBUG_LIBRARY})
    ELSE()
      MESSAGE(WARNING "\n"
        "deal.II was configured with CMAKE_BUILD_TYPE=Debug but no debug tbb\n"
        "library was found. The regular tbb library will be used instead.\n\n"
        )
      LIST(APPEND deal_ii_external_libraries ${TBB_LIBRARY})
    ENDIF()

  ELSE()
    LIST(APPEND deal_ii_external_libraries ${TBB_LIBRARY})
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

  INCLUDE_DIRECTORIES(
    ${CMAKE_SOURCE_DIR}/contrib/tbb/tbb30_104oss/include
    )

  ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/contrib/tbb/tbb30_104oss/src)

  LIST(APPEND deal_ii_additional_object_files
    $<TARGET_OBJECTS:obj_tbb>
    )

ENDMACRO()


CONFIGURE_FEATURE(TBB)

