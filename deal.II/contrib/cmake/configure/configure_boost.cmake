#
# Configuration for the boost library:
#

# Always true. We need it :-]
SET(DEAL_II_WITH_BOOST "ON"
  CACHE STRING "Build deal.II with support for boost." FORCE
  )


MACRO(FIND_FEATURE_BOOST_EXTERNAL var)

  FIND_PACKAGE (Boost COMPONENTS serialization thread)

  IF(Boost_THREAD_FOUND AND Boost_SERIALIZATION_FOUND)
    SET(${var} TRUE)

    # Get rid of this annoying unimportant variable:
    MARK_AS_ADVANCED(Boost_DIR)
  ENDIF()

ENDMACRO()


MACRO(CONFIGURE_FEATURE_BOOST_EXTERNAL var)

  INCLUDE_DIRECTORIES (${Boost_INCLUDE_DIR})

  IF (CMAKE_BUILD_TYPE MATCHES "Debug")
    LIST(APPEND deal_ii_external_libraries
      ${Boost_THREAD_LIBRARY_DEBUG} ${Boost_SERIALIZATION_LIBRARY_DEBUG}
      )
  ELSE()
    LIST(APPEND deal_ii_external_libraries
      ${Boost_THREAD_LIBRARY} ${Boost_SERIALIZATION_LIBRARY}
      )
  ENDIF()

  SET(${var} TRUE)
ENDMACRO()


SET(HAVE_CONTRIB_FEATURE_BOOST TRUE)


MACRO(CONFIGURE_FEATURE_BOOST_CONTRIB var)

  # compile the necessary parts of boost out of ./contrib

  # We need some definitions to use the headers of the bundled boost
  # library:
  ADD_DEFINITIONS("-DBOOST_NO_HASH" "-DBOOST_NO_SLIST")

  INCLUDE_DIRECTORIES(
    ${CMAKE_SOURCE_DIR}/contrib/boost-1.49.0/include
    )

  ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/contrib/boost-1.49.0/libs/serialization/src)

  LIST(APPEND deal_ii_additional_object_files
    $<TARGET_OBJECTS:obj_boost_serialization>
    )

  IF( DEAL_II_USE_MT AND NOT DEAL_II_CAN_USE_CXX1X)
    # If the C++ compiler doesn't completely support the C++1x standard (and
    # consequently we can't use std::thread, std::mutex, etc), then include all
    # the files that form BOOST's thread implementation so that we don't have to
    # build BOOST itself only to get at this small part of it. it also ensures
    # that we use the correct compiler and flags

    ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/contrib/boost-1.49.0/libs/thread/src)

    LIST(APPEND deal_ii_additional_object_files
      $<TARGET_OBJECTS:obj_boost_thread>
      )
  ENDIF()

  SET(${var} TRUE)
ENDMACRO()


MACRO(CONFIGURE_FEATURE_BOOST_ERROR_MESSAGE)
    MESSAGE(SEND_ERROR "
Could not find the boost library!

Please ensure that the boost library is installed on your computer.
If the library is not at a default location, either provide some hints
for the autodetection, or set the relevant variables by hand in ccmake.

Alternatively you may choose to compile the bundled contrib library of
boost by setting DEAL_II_ALLOW_CONTRIB=on or
DEAL_II_FORCE_CONTRIB_BOOST=on.

")
ENDMACRO()


CONFIGURE_FEATURE(BOOST)
