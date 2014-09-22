## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2014 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

#
# Configuration for thread support in deal.II with the help of the tbb
# library:
#

#
# Set up general threading:
# The macro will be included in CONFIGURE_FEATURE_THREADS_EXTERNAL/BUNDLED.
#
MACRO(SETUP_THREADING)
  #
  # Unfortunately the FindThreads macro needs a working C compiler
  #
  IF(CMAKE_C_COMPILER_WORKS)
    #
    # Clear the test flags because FindThreads.cmake will use a C compiler:
    #
    CLEAR_CMAKE_REQUIRED()

    SWITCH_LIBRARY_PREFERENCE()
    FIND_PACKAGE(Threads)
    SWITCH_LIBRARY_PREFERENCE()

    RESET_CMAKE_REQUIRED()

  ELSE()

    #
    # We have no way to query for thread support. Just assume that it is
    # provided by Pthreads...
    #
    MESSAGE(STATUS
      "No suitable C compiler was found! Assuming threading is provided by Pthreads."
      )
    SET_IF_EMPTY(Threads_FOUND TRUE)
    SET_IF_EMPTY(CMAKE_THREAD_LIBS_INIT "-lpthread")
    SET_IF_EMPTY(CMAKE_USE_PTHREADS_INIT TRUE)
  ENDIF()

  IF(NOT Threads_FOUND)
    #
    # TODO: This is a dead end. Threading might be set up with internal TBB
    # so we have no way of returning unsuccessfully...
    #
    MESSAGE(FATAL_ERROR
      "\nInternal configuration error: No Threading support found\n\n"
      )
  ENDIF()

  MARK_AS_ADVANCED(pthread_LIBRARY)

  #
  # Change -lphtread to -pthread for better compatibility on non linux
  # platforms:
  #
  IF("${CMAKE_THREAD_LIBS_INIT}" MATCHES "-lpthread")
    CHECK_CXX_COMPILER_FLAG("-pthread"
      DEAL_II_HAVE_FLAG_pthread
      )
    IF(DEAL_II_HAVE_FLAG_pthread)
      STRING(REPLACE "-lpthread" "-pthread" CMAKE_THREAD_LIBS_INIT
        "${CMAKE_THREAD_LIBS_INIT}"
        )
    ENDIF()
  ENDIF()

  ADD_FLAGS(THREADS_LINKER_FLAGS "${CMAKE_THREAD_LIBS_INIT}")

  #
  # Set up some posix thread specific configuration toggles:
  #
  IF(NOT CMAKE_SYSTEM_NAME MATCHES "Windows")

    IF(NOT CMAKE_USE_PTHREADS_INIT)
      MESSAGE(FATAL_ERROR
        "\nInternal configuration error: Not on Windows but posix thread support unavailable\n\n"
        )
    ENDIF()

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
    RESET_CMAKE_REQUIRED()
    IF(NOT DEAL_II_HAVE_MT_POSIX_BARRIERS)
      SET(DEAL_II_USE_MT_POSIX_NO_BARRIERS TRUE)
    ENDIF()

  ELSE()

    #
    # Poor Windows:
    #
    SET(DEAL_II_USE_MT_POSIX FALSE)
    SET(DEAL_II_USE_MT_POSIX_NO_BARRIERS TRUE)
  ENDIF()

ENDMACRO()


#
# Set up the tbb library:
#

MACRO(FEATURE_THREADS_FIND_EXTERNAL var)
  FIND_PACKAGE(TBB)

  IF(TBB_FOUND)
    SET(${var} TRUE)
  ENDIF()
ENDMACRO()


MACRO(FEATURE_THREADS_CONFIGURE_EXTERNAL)

  IF(CMAKE_BUILD_TYPE MATCHES "Debug")
    IF(TBB_WITH_DEBUG_LIB)
      LIST(APPEND THREADS_DEFINITIONS_DEBUG "TBB_USE_DEBUG" "TBB_DO_ASSERT=1")
      LIST(APPEND THREADS_USER_DEFINITIONS_DEBUG "TBB_USE_DEBUG" "TBB_DO_ASSERT=1")
    ENDIF()
  ENDIF()

  #
  # Workaround for an issue with C++11 mode, non gcc-compilers and missing
  # template<typename T> std::ist_trivially_copyable<T>
  #
  IF( DEAL_II_WITH_CXX11 AND
      NOT DEAL_II_HAVE_CXX11_IS_TRIVIALLY_COPYABLE AND
      NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU" )
    LIST(APPEND THREADS_DEFINITIONS "TBB_IMPLEMENT_CPP0X=1")
    LIST(APPEND THREADS_USER_DEFINITIONS "TBB_IMPLEMENT_CPP0X=1")
  ENDIF()

  SETUP_THREADING()

  LIST(APPEND THREADS_LIBRARIES ${TBB_LIBRARIES})
  LIST(APPEND THREADS_INCLUDE_DIRS ${TBB_INCLUDE_DIRS})
  LIST(APPEND THREADS_USER_INCLUDE_DIRS ${TBB_USER_INCLUDE_DIRS})

ENDMACRO()


MACRO(FEATURE_THREADS_CONFIGURE_BUNDLED)
  #
  # Setup threading (before configuring our build...)
  #
  SETUP_THREADING()

  #
  # We have to disable a bunch of warnings:
  #
  ENABLE_IF_SUPPORTED(THREADS_CXX_FLAGS "-Wno-parentheses")
  ENABLE_IF_SUPPORTED(THREADS_CXX_FLAGS "-Wno-long-long")

  #
  # Add some definitions to use the header files in debug mode:
  #
  IF (CMAKE_BUILD_TYPE MATCHES "Debug")
    LIST(APPEND THREADS_DEFINITIONS_DEBUG "TBB_USE_DEBUG" "TBB_DO_ASSERT=1")
    LIST(APPEND THREADS_USER_DEFINITIONS_DEBUG "TBB_USE_DEBUG" "TBB_DO_ASSERT=1")
  ENDIF()

  #
  # Workaround for an issue with C++11 mode, non gcc-compilers and missing
  # template<typename T> std::ist_trivially_copyable<T>
  #
  IF( DEAL_II_WITH_CXX11 AND
      NOT DEAL_II_HAVE_CXX11_IS_TRIVIALLY_COPYABLE AND
      NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU" )
    LIST(APPEND THREADS_DEFINITIONS "TBB_IMPLEMENT_CPP0X=1")
    LIST(APPEND THREADS_USER_DEFINITIONS "TBB_IMPLEMENT_CPP0X=1")
  ENDIF()

  #
  # tbb uses dlopen/dlclose, so link against libdl.so as well:
  #
  # TODO: Also necessary for external lib, use preference toggle
  #
  LIST(APPEND THREADS_LIBRARIES ${CMAKE_DL_LIBS})

  LIST(APPEND THREADS_BUNDLED_INCLUDE_DIRS ${TBB_FOLDER}/include)
ENDMACRO()


CONFIGURE_FEATURE(THREADS)
